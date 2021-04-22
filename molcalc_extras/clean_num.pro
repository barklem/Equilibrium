function clean_num, s
; cleans out any unwanted stuff around a number

snew = repstr(s, '(', '')
snew = repstr(snew, ')', '')
snew = repstr(snew, '[', '')
snew = repstr(snew, ']', '')
snew = repstr(snew, ' .', '0.')   ; indicates no data
for i = 0, n_elements(snew)-1 do if snew[i] eq '.' then snew[i] = '0.0'
return, snew
end

