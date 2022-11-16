function extract_col, s, j
; extract column j from string s
sout = ' '
if j gt strcnt(s, '<td') then begin
   print, 'EXTRACT_COL WARNING: asks for column which does not exist'
   return, sout
endif
i = 0
cnt = 0
while (cnt ne j) do begin 
   i = strpos(s, '<td', i) 
   if (i ne -1) then begin  
     cnt = cnt + 1 
     i = i + 1  
   endif 
endwhile
ps = strpos(s, '>', i)+1
pe = strpos(s, '</td', ps)
sout = strmid(s, ps, pe-ps)
return, sout
end

function clean_number, s
; cleans out any unwanted stuff around a number

snew = repstr(s, '<sub>', '')
snew = repstr(snew, '</sub>', '')
snew = repstr(snew, '&nbsp;', '.')     ; blank
snew = repstr(snew, 'H<sup>Q</sup>', '')
snew = repstr(snew, 'H', '')
snew = repstr(snew, 'Z', '')
snew = repstr(snew, 'R', '')
snew = repstr(snew, 'V', '')
;snew = repstr(snew, '*', '')
snew = remove_links(snew)
snew = strcompress(snew, /remove_all)
if snew eq '' then snew = '.'
return, snew
end

pro get_symmetry, label, name, Lambda, spin, plusmin, oddeven
; takes the processed label and extracts the symmetries 
; Lambda, 2S+1, +/-, g/u and the state name (e.g. X, a etc)
; Lambda and 2S+1 are returned -1 if not found
; - since +/- is only relevant for Sigma states
;   '-' = -, '+' = +, '.' = both in case of Lambda > 0, or undefined
; - since g/u symmetry is only relevant for homonuclear molecules
;   'g' = g, 'u' = u, '.' = not applicable/undefined

name = ' '
i = strpos(label, '^')
if i ge 0 then begin
   name = strmid(label, 0, i)
   name = repstr(name, '(','')
endif else begin
   i = strpos(label, '(')
   if i ge 0 then begin
      name = strmid(label, 0, i)
   endif else begin
      name = label
   endelse
endelse
name = strtrim(name)

lampos = -1
Lambda = -1
i = strpos(label, 'Sigma')
if i ge 0 then begin
  Lambda = 0
  lampos = i
endif
i = strpos(label, 'Pi')
if i ge 0 then begin
  Lambda = 1
  lampos = i
endif
i = strpos(label, 'Delta')
if i ge 0 then begin
  Lambda = 2
  lampos = i
endif
i = strpos(label, 'Phi')
if i ge 0 then begin
  Lambda = 3
  lampos = i
endif
spin = -1
if lampos gt 0 then begin
  spinpos = strpos(label, '^', lampos, /reverse_search) + 1
  if spinpos gt 0 then begin
    spinstr = strmid(label, spinpos, 2)
    spin = float(spinstr)
  endif  
endif

plusmin = '.'
if Lambda eq 0 then begin
  i = strpos(label, '^+', lampos)
  if i ge 0 then plusmin = '+'
  i = strpos(label, '^-', lampos)
  if i ge 0 then plusmin = '-'
endif

oddeven = '.'
i = strpos(label, '_g', lampos)
if i ge 0 then oddeven = 'g'
i = strpos(label, '_u', lampos)
if i ge 0 then oddeven = 'u'

end




function relabel, s
; converts label string s to something more useful
; needs repstr from astrolib

snew = repstr(s, '&quot;', '"')
snew = repstr(snew, '&#931;', 'Sigma')
snew = repstr(snew, '&#963;', '-sigma')
snew = repstr(snew, '&#928;', 'Pi')
snew = repstr(snew, '&#960;', '-pi')
snew = repstr(snew, '&#916;', 'Delta')
snew = repstr(snew, '&#948;', '-delta')
snew = repstr(snew, '&#934;', 'Phi')
snew = repstr(snew, '&#966;', '-phi')
snew = repstr(snew, '<sup>', '^')
snew = repstr(snew, '</sup>', ' ')
snew = repstr(snew, '<sub>', '_')
snew = repstr(snew, '</sub>', ' ')
snew = repstr(snew, ' )', ')')
snew = remove_links(snew)
snew = strtrim(snew)
return, snew
end

function remove_links, s
snew = s
ps = 0
while ps ge 0 do begin
   ps = strpos(snew, '<a ')
   if ps ge 0 then begin
      pe = strpos(snew, '</a', ps) + 4
      sr = strmid(snew, ps, pe-ps)
      snew = repstr(snew, sr, '')
   endif
endwhile
return, snew
end

function strarrcnt, s, ss
; return number of occurnences of ss in each element of array s
n = n_elements(s)
cnt = intarr(n)*0
for i = 0, n-1 do cnt[i] = strcnt(s[i], ss)
return, cnt
end

 
function strcnt, s, ss
; return number of occurences of ss in s
; s is a single string
i = 0 
cnt = 0 
while (i ne -1) do begin 
   i = strpos(s, ss, i) 
   if (i ne -1) then begin  
     cnt = cnt + 1 
     i = i + 1  
   endif 
endwhile
return, cnt
end


pro parse_mol, infile

; parses html output from NIST
; obtained selecting:
;   - Gas phase ion energetics data
;   - Constants of diatomic molecules
;
; Paul Barklem June 2011

; read html file into string array
nl = file_line(infile)
s = strarr(nl)
openr, lun, infile, /get_lun
readf, lun, s, format ='(a1000)'    ; longest line seen is 600 long
close, lun

; search for ionisation potential
ionpot = 0.
ionpot_err = 0.
ind = strpos(s, 'IE (evaluated)')
indc = where(ind ge 0, nind)
if nind gt 0 then begin
   ; there is often a definition of this statement in the footer, first is data
   s1 = s[indc[0]]                                 ; the string to extract data from
   p1 = ind[indc[0]]                               ; the position of the identified string - data to the right to that
   p2 = strpos(s1, '.', p1)                        ; search forwards for a decimal point
   p3 = strpos(s1, '>', p2, /reverse_search)       ; search backwards for start of entry
   ps = p3 + 1                                     ; start position
   p4 = strpos(s1, ' ', p2)                        ; search forwards for end of entry
   pe = p4 
   ionpot = double(strmid(s1, ps, pe-ps))
 
   ; see if there is an error estimate
   p5 = strpos(s1, '&plusmn;', pe)
   if p5 gt 0 then begin
      ps = strpos(s1, ' ', p5) + 1
      pe = strpos(s1, '<', ps+1)
     ionpot_err = double(strmid(s1, ps, pe-ps))
   endif
endif else begin
   ind = strpos(s, 'IE (eV)')
   indc = where(ind ge 0, nind) + 1
   if nind gt 0 then begin
     s1 = s[indc[0]]
     ipstr = extract_col(s1, 1)
     ; see if there is an error estimate
     p5 = strpos(ipstr, '&plusmn;')
     if p5 gt 0 then begin
       ps = p5 + 8
       ionpot_err = double(strmid(ipstr, ps, 10))
       ionpot = double(strmid(ipstr, 0, ps))
     endif else begin
       ionpot = double(s1)
     endelse
   endif else begin
     print, 'WARNING: Ionisation Energy not found'
   endelse
endelse

print, ionpot, ionpot_err, format = '(2f10.7)'

; now search for molecular states and constants
ind = strpos(s, '<td nowrap="nowrap"')
cnt = strarrcnt(s, '</td>')                 ; counts number of columns
indc = where(ind ge 0 and cnt gt 12, nind)
if nind le 0 then begin
   print, 'STOP: No states found?'
   stop
endif
;print, nind
label = strarr(nind)
Te_str = strarr(nind)
w_str = strarr(nind)
wx_str = strarr(nind)
wy_str = strarr(nind)
B_str = strarr(nind)
alf_str = strarr(nind)
gam_str = strarr(nind)
D_str = strarr(nind)
bet_str = strarr(nind)
re_str = strarr(nind)
for i = 0, nind-1 do begin
    ; extract label in first column
    s1 = s[indc[i]]
    label[i] = relabel(extract_col(s1, 1))
    Te_str[i] = clean_number(extract_col(s1, 2))
    w_str[i] = clean_number(extract_col(s1, 3))
    wx_str[i] = clean_number(extract_col(s1, 4))
    wy_str[i] = clean_number(extract_col(s1, 5))
    B_str[i] = clean_number(extract_col(s1, 6))
    alf_str[i] = clean_number(extract_col(s1, 7))
    gam_str[i] = clean_number(extract_col(s1, 8))
    D_str[i] = clean_number(extract_col(s1, 9))
    bet_str[i] = clean_number(extract_col(s1, 10))
    re_str[i] = clean_number(extract_col(s1, 11))
    get_symmetry, label[i], name, lambda, spin, plusmin, oddeven
    print, label[i], name, lambda, spin, plusmin, oddeven, Te_str[i], w_str[i], wx_str[i], wy_str[i], B_str[i], alf_str[i], gam_str[i], D_str[i], bet_str[i], re_str[i], format = '(a40, a7, i3, i3, a3, a3, 10a12)'
endfor
    


stop
end 
