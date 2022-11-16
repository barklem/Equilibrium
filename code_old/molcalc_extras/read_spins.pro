pro read_spins, names, spins

names = strarr(100)
spins = fltarr(100)
openr, 1, 'nuc_spins.dat'
s = ' '
readf, 1, s, format = '(a30)'   ; comment line
for i = 0, 100 do begin
   s = ' '
   readf, 1, s, format = '(a30)'
   s1 = strsplit(s, ' ', /extract)
   if s1[0] eq 'END' then goto, breakout
   names[i] = s1[0]
   spins[i] = s1[1]
endfor
breakout:
close, 1

names = names[0:i-1]
spins = spins[0:i-1]


end