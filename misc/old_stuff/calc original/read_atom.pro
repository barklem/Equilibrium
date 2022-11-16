pro read_atom, filenm, J, E, ionpot

; reads atom energy levels data from filenm
; returns J and E(cm-1) for each level

 
openr, lun, filenm, /get_lun

; first line is ionisation energy in eV

ionpot = 0.0d0
readf, lun, ionpot

E = dblarr(1000)
j = fltarr(1000)

ns = 0
startp:

  s = ' '
  readf, lun, s, format='(a50)'
  if strmid(s,0,4) eq '    ' then goto, comment  ; nothing in first cols
  if strmid(s,0,1) eq '#' and strmid(s,1,1) ne '#' then goto, comment
  if strmid(s,0,1) eq '-' then goto, comment
  if strmid(s,1,3) eq '---' then goto, comment
  if strmid(s,2,3) eq '---' then goto, comment
  if strmid(s,0,3) eq '###' then goto, endp
  if strpos(s, '-') gt -1 then begin ; a range of J
    s = strsplit(s,' ', /extract)
    s2 = strsplit(s[0],'-', /extract)
    s2b = strsplit(s[2],'[', /extract)  ; remove brackets if any
    s3b = strsplit(s2b,']', /extract)
    s4b = strsplit(s3b[0],'+x', /extract)
    s5b = strsplit(s4b[0],'a', /extract)
    s6b = strsplit(s5b[0],'?', /extract)
    upp = eval_frac(s2[0])
    low = eval_frac(s2[1])
    for i = low, upp, 1 do begin
       e[ns] = double(s6b)
       j[ns] = i
       ns = ns + 1
    endfor
    goto, comment
  endif 
  
  if strpos(s, ',') gt -1 then begin ; a range of J
    s = strsplit(s,' ', /extract)
    s2 = strsplit(s[0],',', /extract)
    s2b = strsplit(s[2],'[', /extract)  ; remove brackets if any
    s3b = strsplit(s2b,']', /extract)
    s4b = strsplit(s3b[0],'+x', /extract)
    s5b = strsplit(s4b[0],'a', /extract)
    s6b = strsplit(s5b[0],'?', /extract)
    upp = eval_frac(s2[0])
    low = eval_frac(s2[1])
    for i = low, upp, 1 do begin
       e[ns] = double(s6b)
       j[ns] = i
       ns = ns + 1
    endfor
    goto, comment
  endif 
  
  
  ; normal case!
  
  
    s = strsplit(s,' ', /extract)
;    print, s[0], '***', s[1], '***', s[2]
    s2 = strsplit(s[2],'[', /extract)  ; remove brackets if any
    s3 = strsplit(s2,']', /extract)
    s4 = strsplit(s3[0],'+x', /extract)
    s5 = strsplit(s4[0],'a', /extract)
    s6 = strsplit(s5[0],'?', /extract)
    if s6 eq '_' or s6 eq '|' then begin
      e[ns] = e[ns-1]
    endif else begin
;      print, '***',s6
      e[ns] = double(s6[0])
    endelse
    j[ns] = eval_frac(s[0])
    ns = ns + 1
    
;  print, j[ns], e[ns]
  comment:
  goto, startp
  endp:
  
free_lun, lun

j = j[0:ns-1]
e = e[0:ns-1]

return
end
