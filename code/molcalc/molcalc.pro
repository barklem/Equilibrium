; The main program

pro molcalc

; turn on this switch to calculate negative molecular ion equilibrium constants 
; for 3 body dissociation i.e. A + B + e <-> AB- , so K = p(A)p(B)p(e)/p(AB-)
; rather than default 2 body i.e. A + B- <-> AB- , so K = p(A)p(B-)/p(AB-)
; note this means data with mixed dimensions and units!!!   

negative3 = 0

; read all data

read_diss, molname, DE, components=molcomp
nmol = n_elements(molname)
molecule_components = strarr(nmol,3)
molecule_components[*,0:1] = molcomp

read_weights, wname, wweight

; file locations 
;atomfile = '../remo/molcalc_test/allatomdata_new.sav'
atomfile = '../remo/molcalc_test/alldata_atomsonly.sav'
atomfile = '../remo/allatomdata_new.sav'
atomfile = '../molcalc_atoms/allatomdata_vlargeTgrid.sav'
atomfile = '../molcalc_atoms/allatomdata_exlargeTgrid.sav'
atomfile = '../molcalc_atoms/allatomdata_post_fix_2022_origT.sav'
atomfile = '../molcalc_atoms/allatomdata_post_fix_2022_exlargeT.sav'
moldir = './Final_Mols_idl/'  ; should have ending slash
outdir1 = './PartFuncs_Results/'
 
; calculate thermodynamic properties on the following T grid

; Temperature grid for Barklem & Collet (2016) paper
T = [1e-5,  1e-4,  1e-3,  1e-2,  0.1,   0.15,  0.2,   0.3,   0.5,   0.7, $
     1.0,   1.3,   1.7,   2.0,   3.0,   5.0,   7.,    10.,   15.,   20., $
     30.,   50.,   70.,   100.,  130.,  170.,  200.,  250.,  300.,  500., $
     700.,  1000., 1500., 2000., 3000., 4000., 5000., 6000., 7000., 8000., $
     9000., 10000. ]

; Temperature grid, Oct. 2016 calculations  ("extended T grid")
;  T = 10.^(findgen(1001)/100.-6)

; Temperature grid, Aug. 2020 and later calculations, 10^-6 to 10^6, 100 / 400 points per decade
;  T = 10.^(findgen(12*100+1)/100.-6)        ;   ("vlarge T grid")
  T = 10.^(findgen(12*400+1)/400.-6)        ;   ("extremely large T grid")

theta = 5040.d0/T
nT = n_elements(T)

; below we perform all tasks on each molecule in turn

Qmol = dblarr(nmol, nT)
Qmol4 = dblarr(nmol, nT)
QmolST = dblarr(nmol, nT)
QmolHH = dblarr(nmol, nT)
lgKpmol = dblarr(nmol, nT)
Dmol = dblarr(nmol)
EZEROmol = dblarr(nmol)

; prepare clean file for output of problem cases
logfile = 'logfile.txt'
openw, lun, logfile, /get_lun
printf, lun, 'Log file for molcalc run'
printf, lun, 'starting ', systime()
close, lun
free_lun, lun

print, ' '
print, 'Calculating for molecules:'
for i = 0, nmol-1 do begin 
;for i = 0, 40 do begin  
  restore, moldir + molname[i] + '_hh.idl'  
  moldat_hh = moldat
  restore, moldir + molname[i] + '.idl'
  Dmol[i] = moldat.D
  
  print, ' '
  print,  string(i+1, '(i3)')+'/'+string(nmol, '(i3)')+': ', molname[i], moldat.ns, moldat.D
  print, 'Partition function: ', molname[i]
  print, 'Calculating with merged data'
  comp_mol_partf, T, moldat, Qi, Qi4=Qi4, QiST=QiST, errorfile=logfile, /showinfo, Ezero=Ezero
  EZEROmol[i] = Ezero 
  print, ' Ezero = ', Ezero 
  print, 'Calculating with HH data'
  comp_mol_partf, T, moldat_hh, Qi_hh, Qi4=Qi4_hh, QiST=QiST_hh
  Qmol[i, *] = Qi
  Qmol4[i, *] = Qi4
  QmolST[i, *] = QiST_hh
  QmolHH[i, *] = Qi_hh
  
   lgKp = dblarr(nT)
   print, 'Equilibrium constant: ', molname[i]
;  i1 = where(atomid eq molatid[i,0], n1) 
;  i2 = where(atomid eq molatid[i,1], n2) 
   
; atomic weights
   atn1 = molcomp[i,0]
   atn1 = repstr(atn1, '+', '')
   atn1 = repstr(atn1, '-', '')
   atn2 = molcomp[i,1]
   atn2 = repstr(atn2, '+', '')
   atn2 = repstr(atn2, '-', '')
   ind = where(wname eq atn1)
   atm1 = wweight(ind)
   ind = where(wname eq atn2)
   atm2 = wweight(ind)
   atm1 = atm1[0]
   atm2 = atm2[0]

; construct species names
   atn_neg = ' '

   atn1_ineg = 0
   atn1 = molcomp[i,0]
   i1 = strpos(atn1, '-')
   i2 = strpos(atn1, '+')
   if (i1 eq -1 and i2 eq -1) then atn1 = atn1 + '_I'
   atn1 = repstr(atn1, '+', '_II')
   if i1 gt 0 then atn1_ineg = 1

   atn2_ineg = 0
   atn2 = molcomp[i,1]
   i1 = strpos(atn2, '-')
   i2 = strpos(atn2, '+')
   if (i1 eq -1 and i2 eq -1) then atn2 = atn2 + '_I'
   atn2 = repstr(atn2, '+', '_II')
   if i1 gt 0 then atn2_ineg = 1

   atn_ineg = 0
   if (atn1_ineg or atn2_ineg) then atn_ineg = 1

   ; if negative ion and 3 body calculation, construct correct data
   if negative3 and atn_ineg then begin
      if atn1_ineg then begin 
         atn_neg = atn1
         atn1 = repstr(atn1, '-', '_I')
         atn1_clean = repstr(atn1, '_I', ' ')
         atn1_clean = strcompress(atn1_clean, /remove_all)
         molecule_components[i,0] = atn1_clean
         molecule_components[i,2] = 'e'
      endif   
      if atn2_ineg then begin 
         atn_neg = atn2
         atn2 = repstr(atn2, '-', '_I')
         atn2_clean = repstr(atn2, '_I', ' ')
         atn2_clean = strcompress(atn2_clean, /remove_all)
         molecule_components[i,1] = atn2_clean
         molecule_components[i,2] = 'e'
      endif   
   endif

   ; make a version with file ending and without
   atn1c = atn1
   atn1 = atn1 + '.sav'
   atn2c = atn2
   atn2 = atn2 + '.sav'
   atn_negc = atn_neg
   atn_neg = atn_neg + '.sav'
   
; use a single collected file  
   
  Tstore = T 
  restore, atomfile
  ind1 = where(atn1c eq atomid)
  ind2 = where(atn2c eq atomid)  
   if ((ind1 ge 0) and (ind2 ge 0)) then begin
      print, 'using ' + atn1c + ' and ' + atn2c
      print, 'masses: ', atm1, atm2
   endif else begin
      if ind1 lt 0 then print, atn1c + ' not found'
      if ind2 lt 0 then print, atn2c + ' not found'
      stop
      goto, skipKp
   endelse
   
   Tatom = T
   Qa1 = reform(Qatom[ind1,*])
   Qa2 = reform(Qatom[ind2,*])
   T = Tstore

   ; if needed extract negative ion ionisation potential for 3 body
   if negative3 and atn_ineg then begin
      ind3 = where(atn_negc eq atomid)
      if (ind3 ge 0) then begin
         print, 'using ' + atn_negc + ' for negative ion'
      endif else begin
        print, atn_negc + ' not found'
        stop
      endelse

      negative_ip = atom_potion[ind3]
      negative_ip = negative_ip[0]
   endif   

; end read atom data
  
  if not array_equal(Tatom, T) then begin
     print, ' WARNING: Atom and molecule temerature grids not the same, interpolating and maybe extrapolating'
     read, ' continue?  (y/other)', ans
     if ans ne 'y' then stop
  endif 

  ; put on same T grid
  Qa1 = interpol(Qa1, Tatom, T)
  Qa2 = interpol(Qa2, Tatom, T)
  
  mu = (atm1 * atm2) / (atm1 + atm2)

  if negative3 and atn_ineg then begin
     comp_mol_equil3, T, moldat, Qa1, Qa2, Qi, mu, negative_ip, lgKp
     Dmol[i] = Dmol[i] + negative_ip
  endif else begin 
     comp_mol_equil, T, moldat, Qa1, Qa2, Qi, mu, lgKp
  endelse   
   
  lgKpmol[i, *] = lgKp
  
  mol_components = molecule_components[i,*]
  save, file = outdir1 + molname[i]+'.idl', lgKp, Qi, T, mol_components   ; save partition func and equil constant

;stop

skipKp:
  
; save all data in a single file

;print, ' '
;print, 'saving all data'
;molid = molname
;save, file = outdir1 + 'allpartf.idl', T, Qmol, Qmol4, QmolSt, QmolHH, molid, Dmol

endfor

; this is file with all atomic data with molecular parts cut out by hand
Tstore = T
;restore, '../remo/molcalc_test/allatomdata_new.sav'
restore, atomfile
Tatom = T
T = Tstore
num_atm = n_elements(atomid)
Qatom_new = dblarr(num_atm, nT)
for i = 0, num_atm-1 do Qatom_new[i,*] = interpol(Qatom[i,*], Tatom, T)
Qatom = Qatom_new



openw, lun, logfile, /get_lun, /append
printf, lun, systime()
close, lun
free_lun, lun

; save all data in a single file

print, ' '
print, 'saving all data'
molid = molname
save, file = outdir1 + 'allatom.idl', T, Qatom, atomid, atom_potion
save, file = outdir1 + 'allpartf.idl', T, Qmol, Qmol4, QmolSt, QmolHH, molid, Dmol, EZEROmol
save, file = outdir1 + 'allequil.idl', T, lgKpmol, molid, Dmol, molecule_components

print, ' '
print, 'normal end'
print, ' '
stop
end


