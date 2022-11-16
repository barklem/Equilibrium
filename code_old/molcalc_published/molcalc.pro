; The main program
;
; Reads in list of atoms and molecules to be calculated 
; with files for fundamental data
; stored in designated file mollist
;
; Note, the program is kept general to allow molecules larger than 
; diatomic, but they are not yet supported in the calculation routines
; 
; Paul Barklem, Feb 2007
; Major revision, Nov 2011


pro molcalc;, mollist

; read all data

read_diss, molname, DE, components=molcomp
nmol = n_elements(molname)

read_weights, wname, wweight


; file locations 
atomdir = '../remo/molcalc_test/atompf/' 
;atomfile = '../remo/molcalc_test/allatomdata_new.sav'
atomfile = '../remo/molcalc_test/alldata_atomsonly.sav'
moldir = './Final_Mols_idl/'  ; should have ending slash
outdir1 = './PartFuncs_Results/'
 
; calculate thermodynamic properties on the following T grid

;T = 10.^(findgen(101) * 0.06 - 2.)
;T = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7., 10., 20., 30., 50., 70., 100., 200., 300., 500., 700., 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000. ]
; this is a temperature grid built by comparison with Jeff Valenti's adaptive grid based on the fine grid results
; which is built to allow cubic spline with accuracy better than 1e-4.  i.e. we build a "even-numbered" grid similar to his
T = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.7, 2.0, 3.0, 5.0, 7., 10., 15., 20., 30., 50., 70., $
     100., 130., 170., 200., 250., 300., 500., 700., 1000., 1500., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000. ]
;T = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (findgen(10000)+1.)/10000. * 10000. ]
;T = 10.^(findgen(501)/100-1.)
;T = 10.^(findgen(1001)/100.-5)

;minT = 1e-5
;maxT = 1.e4
;nTT = 301   ; choose a number like 101 to hit round numbers best
;lgT = alog10(minT) + findgen(nTT) / (nTT - 1) * (alog10(maxT) - alog10(minT))
;T = 10.^lgT

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
   atn1 = molcomp[i,0]
   i1 = strpos(atn1, '-')
   i2 = strpos(atn1, '+')
   if (i1 eq -1 and i2 eq -1) then atn1 = atn1 + '_I'
   atn1 = repstr(atn1, '+', '_II')
   atn1c = atn1
   atn1 = atn1 + '.sav'
   atn2 = molcomp[i,1]
   i1 = strpos(atn2, '-')
   i2 = strpos(atn2, '+')
   if (i1 eq -1 and i2 eq -1) then atn2 = atn2 + '_I'
   atn2 = repstr(atn2, '+', '_II')
   atn2c = atn2
   atn2 = atn2 + '.sav'
   
iatomfiles = 0   ; use individual files or use a collected files   
   
if iatomfiles then begin   
   q1 = file_test(atomdir + atn1)
   q2 = file_test(atomdir + atn2)
   if (q1 and q2) then begin
      print, 'using ' + atn1 + ' and ' + atn2
      print, 'masses: ', atm1, atm2
   endif else begin
      if q1 eq 0 then begin
         print, atn1 + ' not found'
         printf, lun, atn1 + ' not found - skipping ', molname[i]
      endif
      if q2 eq 0 then begin
         print, atn2 + ' not found'
         printf, lun, atn2 + ' not found - skipping ', molname[i]
      endif
      
      goto, skipKp
   endelse
  
  Qi_store = Qi
  Tstore = T 
  restore, atomdir + atn1  
  Qa1 = QI
  restore, atomdir + atn2 
  Qa2 = QI
  Tatom = T
  T = Tstore
  Qi = Qi_store
  
endif else begin  
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
   Qa1 = Qatom[ind1,*]
   Qa2 = Qatom[ind2,*]
   T = Tstore
endelse
  
  ; put on same T grid
  Qa1 = interpol(Qa1, Tatom, T)
  Qa2 = interpol(Qa2, Tatom, T)
  
  mu = (atm1 * atm2) / (atm1 + atm2)
  comp_mol_equil, T, moldat, Qa1, Qa2, Qi, mu, lgKp
   lgKpmol[i, *] = lgKp
  
  save, file = outdir1 + molname[i]+'.idl', lgKp, Qi, T   ; save partition func and equil constant

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
save, file = outdir1 + 'allequil.idl', T, lgKpmol, molid, Dmol

print, ' '
print, 'normal end'
print, ' '
stop
end


