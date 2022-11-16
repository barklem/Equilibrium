restore, 'allatomdata_post_fix_2022_origT.sav'
Qatom_new = Qatom
restore, 'allatomdata_pre_fix_2022_origT.sav'

nT = n_elements(T)
nA = n_elements(atomid)

iT = nT-1  ; 10000 K
iT = 36    ; 5000 K

for i = 0, nA-1 do begin
	diff = (Qatom_new[i, iT] - Qatom[i, iT]) / Qatom[i, iT]
	if abs(diff) gt 0.d0 then print, atomid[i], diff, format = '(a5, e12.2)'
endfor
stop
end	
