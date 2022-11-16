function eval_frac, s

; evaluates a fraction string like '5/2'
; or changes any number string to float

if strpos(s, '/') gt -1 then begin ; have a frac
  s = strsplit(s,'/', /extract)
  num = float(s[0])
  den = float(s[1])
  return, num/den
endif else begin
  return, float(s)
endelse

end

