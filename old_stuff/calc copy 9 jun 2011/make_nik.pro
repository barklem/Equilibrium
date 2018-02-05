; adjusts alldata so negative ions have negative ionisation potentials
; to fit Niks format

restore, 'alldata.sav'

ind = where(strpos(atomid, '-') gt 0, nind)
atom_potion[ind] = -atom_potion[ind]

save, file = 'alldata_neg.sav', T, lgKpmol, Qmol, Qatom, molid, atomid, atom_potion, Dmol

end
