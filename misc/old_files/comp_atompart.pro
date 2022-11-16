pro comp_atompart, filenm, T, Qi

; computes partition function for atom
; data from filenm
; internal partition func Qi returned for temps T

 
openr, 1, filenm
Ilim = 0.
s = ' '
readf, 1, s ; comment line
readf, 1, Ilim
readf, 1, s ; comment line

E = dblarr(1000)
j = fltarr(1000)

ns = 0
startp:

  s = ' '
  readf, 1, s, format='(a50)'
  if strmid(s,0,4) eq '    ' then goto, comment  ; nothing in first cols
  if strmid(s,0,1) eq '#' and strmid(s,1,1) ne '#' then goto, comment
  if strmid(s,0,1) eq '-' then goto, comment
  if strmid(s,0,3) eq '###' then goto, endp
  if strpos(s, '-') gt -1 then begin ; a range of J
    s = strsplit(s,' ', /extract)
    s2 = strsplit(s[0],'-', /extract)
    s2b = strsplit(s[2],'[', /extract)  ; remove brackets if any
    s3b = strsplit(s2b,']', /extract)
    s4b = strsplit(s3b[0],'+x', /extract)
    upp = eval_frac(s2[0])
    low = eval_frac(s2[1])
    for i = low, upp, 1 do begin
       e[ns] = double(s4b)
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
    upp = eval_frac(s2[0])
    low = eval_frac(s2[1])
    for i = low, upp, 1 do begin
       e[ns] = double(s4b)
       j[ns] = i
       ns = ns + 1
    endfor
    goto, comment
  endif 
  
  
  ; normal case!
  
  
    s = strsplit(s,' ', /extract)
    s2 = strsplit(s[2],'[', /extract)  ; remove brackets if any
    s3 = strsplit(s2,']', /extract)
    s4 = strsplit(s3[0],'+x', /extract)
    if s4 eq '_' then begin
      e[ns] = e[ns-1]
    endif else begin
      e[ns] = double(s3[0])
    endelse
    j[ns] = eval_frac(s[0])
    ns = ns + 1
    
;  print, j[ns], e[ns]
  comment:
  goto, startp
  endp:
  
close, 1

g = 2*j+1
Ej = E /8065.54 * 1.602d-19  ; cm-1 > J
k = 1.38066d-23

Qi = 0.d0 * T

for i = 0, ns -1 do begin
  Qi = Qi + g[i]*exp(-Ej[i]/k/T)
;  print, g[i], exp(-Ej[i]/k/T), Qi[0]
endfor

return
end
