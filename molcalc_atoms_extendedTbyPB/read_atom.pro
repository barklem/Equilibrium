pro read_atom, fname, J, E, ionpot, debug=debug, nlev=nlev, ns=ns

;+
; NAME: READ_ATOM
; PURPOSE: read atomic energy levels data from ASCII file
; INPUT: fname = file name
; OUTPUT: Energies E (cm^-1) and angular momenta J for all tabulated levels
;-
; Bug fix on reading J values where more than 1, see line 170 - PB
; 
  IdStr = 'READ_ATOM'
  
; default keywords
  if keyword_set(debug) then do_debug = 1 else do_debug = 0
  if n_elements(nlev) eq 0 then nlev=4000L

; initialise variables and keywords
  is_spin_half = 0
  is_end_str = 0
  is_comment = 0
  is_header = 0
  E  = dblarr(nlev)             ; energy (cm^-1)
  j  = fltarr(nlev)             ; angular momentum, quantum number
  term = lonarr(nlev)           ; term number
  iterm = -1                    ; term counter (-1: no terms, 0: first term, etc.)
  ionstage = lonarr(nlev)       ; ionisation stage
  iion = -1                     ; ionisation stage counter
  ns = 0                        ; actual number of levels, counter
  ionpot = 0.                   ; ionisation energy (eV !)

; format input string
  fmt_input_str = '(A127)'

  
; open file and read ionisation energy from first line
  openr, lun, /get_lun, fname
  readf, lun, ionpot

; loop over lines of text in file
  while ~EOF(lun) and ~is_end_str do begin

; reset input string s
; read current line and store it in s
; copy string s to sbak
     s = ' '
     readf, lun, s, format=fmt_input_str
     sbak = s
     if (do_debug) then print, '% '+IdStr+': read string', s
     
; check if header string 
     is_header_j = strpos(s,'J') gt -1
     is_header_e = strpos(s,'E') gt -1
     is_header = is_header_j and is_header_e
     if (is_header and do_debug) then print, '% '+IdStr+': is_header_j, is_header_e=',is_header_j,is_header_e

     if ~is_header then begin
        
; check if comment string
        is_comment = (strmid(s,0,1) eq '#' and strmid(s,1,1) ne '#') or strmid(s,0,1) eq '-'
        if (is_comment and do_debug) then  print, '% '+IdStr+': comment string'

; check if end string
        is_end_str = strmid(s,0,3) eq '###'
        if (is_end_str and do_debug) then  print, '% '+IdStr+': end string'

        
; if not comment nor end string, proceed with string processing:
; all other strings are expected to contain two fields
; delimited by '|'
        
        if ~(is_comment or is_end_str) then begin

; split s string in two, J part and E part
           res   = strsplit(s,'|',/extract)
           str_j = res[0]
           str_e = res[1]
           
; remove blanks 
           str_j = strcompress(str_j, /remove_all)
           str_e = strcompress(str_e, /remove_all)

; if processed str_j and str_e strings are empty, then this marks
; a term separator: increase term counter
           is_new_term = (str_j eq '' and str_e eq '')
           if is_new_term then iterm = iterm + 1
           
; if there is at least one term, and this is the first term occurring,
; then set iion to 0
           if is_new_term and iion eq -1 then iion = 0
           
; check if strings only contain '-' characters
; it could either be a simple separator or mark the ionisation energy (---);
; either way, we ignore those entries for now and treat any levels
; beyond that as autoionisation levels

           is_ionpot_str = str_j eq '---'
           is_separator = strjoin(strsplit(str_j,'-',/extract)) eq '' and $
                          strjoin(strsplit(str_e,'-',/extract)) eq ''

; if ionpot string, increase iion counter
           if is_ionpot_str then iion = iion + 1

; if not a ionisation potential string or separator marker either,
; then treat it as atomic energy level data
           
           if ~(is_ionpot_str or is_separator) then begin
              
; remove any left-over characters from HTML trimming of NIST output
              str_j_tmp = strsplit(str_j,'[()]+ax?', /extract) ; remove brackets if any
              str_e_tmp = strsplit(str_e,'[()]+ax?', /extract) ; remove brackets if any
              str_j     = strjoin(str_j_tmp)
              str_e     = strjoin(str_e_tmp)

              
; is there a range of J values?
              
; case 1: list of J values separated by or
              has_or = (stregex(str_j,'or') gt -1)
              if (has_or) then begin
                 if (do_debug) then print, '% '+IdStr+': string has a range of Js'
                 str_j = strsplit(str_j,'or',/regex,/extract,count=nsplit)
                 jtmp  = fltarr(nsplit)    
                 for i=0L,nsplit-1 do begin
                    jtmp[i] = eval_frac(str_j[i])
                 endfor
; if energy information is missing, then mark it with -1.0 for now
                 if str_e eq '' or str_e eq '_' then begin
                    print, '% '+IdStr+': Warning! Missing energy information!'
                    E[ns] = -1.0
                 endif else begin
                    E[ns] = double(str_e)
                 endelse
                 jmean = mean(jtmp)
                 if (do_debug) then $
                    print, '% '+IdStr+': Range of Js: '+str_j+', Jmean =', jmean
                 j[ns] = jmean
                 term[ns] = iterm
                 ionstage[ns] = iion
                 ns = ns + 1
              endif 

;; check if str_j contains only '---': 
;; to be on the safe side, assign J=1/2 for half-integer spin case and
;; J=0 for integer case 
;  if (str_j eq '---') then begin
;     if (is_half_spin) then j[ns]=0.5 else j[ns]=0.0
;     e[ns] = double(str_e)
;     ns = ns + 1
;     goto, comment
;  endif
              
; case 2: range of J values, separated by commas or dashes
; note: case str_j='---' is excluded by if statement above;
; note: there appears not to be situations where energy
;       information is missing in the data for this particular case

              if (~has_or) then begin
                 has_dash  = (strpos(str_j, '-') gt -1)
                 has_comma = (strpos(str_j, ',') gt -1)
                 if (has_dash or has_comma) then begin
                    if (do_debug) then $
                       print, '% '+IdStr+': string has a range of Js separated by commas or dashes'
; first, replace all dashes with commas   PB2022 - I don't find any such cases
                    str_j  = strjoin(strsplit(str_j, '-', /extract), ',')
; then separate substrings assuming commas as separator
                    str_j_arr = strsplit(str_j, ',', /extract, count=nsplit)
                    for i=0L,nsplit-1 do begin
                       E[ns] = double(str_e)
                       j[ns] = double(str_j_arr[i])   ;BUG FIX PB2022 - missing eval_frac function
                       ;j[ns] = double(eval_frac(str_j_arr[i]))
                       term[ns] = iterm
                       ionstage[ns] = iion
                       ns = ns + 1
                    endfor
                 endif 
              endif             ; ~has_or
; finally, normal case
              if ~(has_or or has_dash or has_comma) then begin
                 if (do_debug) then $
                    print, '% '+IdStr+': normal case'
; if energy information is missing, then mark it with -1.0 for now
                 if str_e eq '' or str_e eq '_' then begin
                    print, '% '+IdStr+': Warning! Missing energy information!'
                    E[ns] = -1.0
                 endif else begin
                    E[ns] = double(str_e[0])
                 endelse
                 j[ns] = eval_frac(str_j)
                 term[ns] = iterm
                 ionstage[ns] = iion
                 ns = ns + 1
              endif             ; normal case
              
           endif                ; ~(is_ionpot_str or is_separator)
        endif                   ; ~(is_comment or is_end_str)
     endif                      ; ~is_header

; check if spin is integer or half integer
     if (ns eq 1) then begin
        if round(j[0]) eq j[0] then is_half_spin=0 else is_half_spin=1
     endif

  endwhile
  
; close file
  free_lun, lun

; print number of lines
  if (do_debug) then print, '% '+IdStr+': number of levels: ns=', ns

; reform arrays and exit
  J = reform(j[0:ns-1])
  E = reform(E[0:ns-1])
  term = reform(term[0:ns-1])
  ionstage = reform(term[0:ns-1])

; take a copy of E array
  Ebak = E
  
; deal with missing energy entries (E=-1)
  ww_neg = where(Ebak eq -1., nww_neg)

  if (nww_neg gt 0) then begin
; loop over missing entries
     for i=0,nww_neg-1 do begin
        k = ww_neg[i]
; find valid entries belonging to the same term
        ww = where(term eq term[k] and Ebak ne -1.0, nww)

; normal case: more than two valid entries in the same term
; sort entries in order of increasing J;
; interpolate (or extrapolate) E linearly in J
        if (nww ge 2) then begin
           wws = sort(j[ww])
           E[k] = interpol(Ebak[ww[wws]],J[ww[wws]],J[k])
        endif
; case nww eq 1:        
; use same E as valid entry
        if (nww eq 1) then begin
           E[k] = Ebak[ww[0]]
        endif
; case nww eq 0:
; in lack of a better option, use previous value in E array
; if very first entry, set E=0
        if (nww eq 0) then begin
           if (k gt 0) then E[k] = E[k-1] else E[k]=0.
        endif
     endfor
  endif

;  stop
  
  return
end
