pro read_weights, names, weights

names = strarr(130)
weights = fltarr(130)
openr, 1, 'atomic_weights.dat'
s = ' '
;readf, 1, s, format = '(a50)'   ; comment line
for i = 0, 130 do begin
   s = ' '
   readf, 1, s, format = '(a50)'
   s1 = strsplit(s, ' ', /extract)
   if s1[0] eq 'END' then goto, breakout
   names[i] = s1[1]
   weights[i] = s1[3]
endfor
breakout:
close, 1

names = strcompress(names[0:i-1], /remove_all)
weights = weights[0:i-1]


end