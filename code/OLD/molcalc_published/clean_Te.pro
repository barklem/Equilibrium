function clean_Te, s
; cleans out any unwanted stuff around a number

snew = repstr(s, '(', '')
snew = repstr(snew, ')', '')
snew = repstr(snew, '[', '')
snew = repstr(snew, ']', '')
; handle case of Te in eV
for i = 0, n_elements(snew)-1 do begin
  iev = strpos(snew[i], 'eV')
  if iev ge 0 then begin
     snew[i] = repstr(snew[i], 'eV', '')
     snew[i] = string(double(snew[i])*8065.54, '(f10.1)')
  endif
endfor


snew = repstr(snew, ' .', '-1.0')   ; indicates no data
for i = 0, n_elements(snew)-1 do if snew[i] eq '.' then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if snew[i] eq '' then snew[i] = '-1.0'

; must deal with cases like x1, and x1+12345., i.e. relative energies

for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'x') ne -1L then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'y') ne -1L then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'z') ne -1L then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'a') ne -1L then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'b') ne -1L then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'c') ne -1L then snew[i] = '-1.0'
for i = 0, n_elements(snew)-1 do if strpos(snew[i], 'd') ne -1L then snew[i] = '-1.0'

return, snew
end