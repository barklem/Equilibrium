pro read_ionpot_nist, inputfile, atomnum_arr, elemsym_arr, ionsym_arr, ionpot_arr, $
                      dir=dir, maxitems=maxitems, verbose=verbose, $
                      nitems=nitems, is_theor_arr=is_theor_arr, is_interp_arr=is_interp_arr
;+
;  NAME: READ_IONPOT_NIST
;
;  PURPOSE: read file with ionisation potentials and sort data into
;  various arrays
;  OPTIONS (INPUT):
;    dir      = directory where input file is stored
;    maxitems = max allowed size of output arrays
;    verbose  = verbose output?
;-

; some default values  
  if n_elements(dir) eq 0 then dir = '.'
  if n_elements(maxitems) eq 0 then maxitems = 500L
  if n_elements(verbose) eq 0 then verbose = 0L

; initialise temporary arrays
  tmpAtomNumArr = lonarr(maxitems)
  tmpElemSymArr = strarr(maxitems)
  tmpIonSymArr  = strarr(maxitems)
  tmpIonPotArr  = fltarr(maxitems)
  tmpIsTheorArr = lonarr(maxitems)
  tmpIsInterpArr= lonarr(maxitems)
  
; initialise number of items
  nitems = 0L

; initialise counter
  iitem = 0L
  
; main program starts here
  cd, current=current_dir

; open file
  cd, dir
  openr,lun,/get,inputfile

; loop over lines until end-of-file
  while ~EOF(lun) do begin

; initialise string for input
; initialise logical switches
     s = ' '
     is_comment = 0
     is_header  = 0
     is_data    = 0
; initialise output variable (format)
     atomicNum  = 0
     elemSym    = ''
     ionSym     = ''
     ionPot     = 0.0
; logical switches: theoretical or interpolated value?
     is_theor   = 0
     is_interp  = 0
; left and right bracket symbols, for printout
     left_bracket  = ''
     right_bracket =''
     
; read string s
     readf,lun,s
     
; check if string contains data or else (comment or header)
     if strmid(s,0,1) eq '-' then is_comment = 1 else is_comment = 0
     if ~is_comment then begin
        if strmid(s,0,3) eq 'At.' then is_header = 1 else is_header = 0
     endif
     if ~(is_comment or is_header) then is_data = 1 else is_data = 0

; if data, then process string to extract information
     if is_data then begin

; split into substrings, '|' separator
        subStrArr = strsplit(s,'|',/extract)

; extract (read) atomic number from first substring
        reads, subStrArr[0], atomicNum

; extract element and ion symbols from second string
        subStrArr1 = strsplit(subStrArr[1],' ',/extract)
        elemSym = subStrArr1[0]
        ionSym  = subStrArr1[1]

; check if theoretical or interpolated value
        if strpos(subStrArr[2],'(') gt -1 then begin
           is_theor = 1
           left_bracket  = '('
           right_bracket = ')'
        endif
        if strpos(subStrArr[2],'[') gt -1 then begin
           is_interp = 1
           left_bracket  = '['
           right_bracket = ']'
        endif
; remove brackets, if any:
        subStr2 = strsplit(subStrArr[2],'[()]+ax?', /extract)
; extract ionisation energy
        reads, subStr2, ionpot

; store results into arrays
        tmpAtomNumArr[iitem] = atomicNum
        tmpElemSymArr[iitem] = elemSym
        tmpIonSymArr[iitem]  = ionSym
        tmpIonPotArr[iitem]  = ionpot
        tmpIsTheorArr[iitem] = is_theor
        tmpIsInterpArr[iitem]= is_interp
; print results?
        if (verbose) then print, atomicNum,':',elemSym,ionSym,left_bracket,ionPot,right_bracket,format='(2X,I3,X,A1,2X,A2,X,A3,3X,A1,F7.4,X,A1)'

;increase counter
        iitem = iitem + 1
     endif

; total number of valid items
     nitems = iitem
  endwhile

  if nitems gt 0 then begin
     atomnum_arr   = reform(tmpAtomNumArr[0:nitems-1])
     elemsym_arr   = reform(tmpElemSymArr[0:nitems-1])
     ionsym_arr    = reform(tmpIonSymArr[0:nitems-1])
     ionpot_arr    = reform(tmpIonPotArr[0:nitems-1])
     is_theor_arr  = reform(tmpIsTheorArr[0:nitems-1]) 
     is_interp_arr = reform(tmpIsInterpArr[0:nitems-1])
  endif else begin
     atomnum_arr   = 0L
     elemsym_arr   = ''
     ionsym_arr    = ''
     ionpot_arr    = -1.0
     is_theor_arr  = 0L
     is_interp_arr = 0L
  endelse

; close file
  free_lun,lun
; end
  cd, current_dir
end
