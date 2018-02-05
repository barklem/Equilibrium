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

pro molcalc, mollist

; read all data

openr, lun, mollist, /get_lun

; read atom list first

maxnatom = 100
atomid = strarr(maxnatom)
atomdata = strarr(maxnatom)

natom = 0

for i = 0, maxnatom-1 do begin
  s1 = ' '
  readf, lun, s1
  s1 = strcompress(s1)
  s = strsplit(s1, ' ', /extract)
  
  if (strmid(s[0],0,1) eq '#') then goto, skipatom
  if (strmid(s[0],0,3) eq 'end') then goto, endatoms
  atomid[natom] = s[0]
  atomdata[natom] = s[1]
  natom = natom + 1
  skipatom:
endfor

endatoms:
atomid = atomid[0:natom-1]
atomdata = atomdata[0:natom-1]
 
; read molecule list now

maxnmol = 100
maxnmolat = 2
molid = strarr(maxnmol)
moldata = strarr(maxnmol)
molnat = intarr(maxnmol)
molatid = strarr(maxnmol, maxnmolat)
molatmass = fltarr(maxnmol, maxnmolat)

nmol = 0

for i = 0, maxnmol-1 do begin
  s1 = ' '
  readf, lun, s1
  s1 = strcompress(s1)
  s = strsplit(s1, ' ', /extract)
  
  if (strmid(s[0],0,1) eq '#') then goto, skipmol
  if (strmid(s[0],0,3) eq 'end') then goto, endmolecules
  molid[nmol] = s[0]
  moldata[nmol] = s[1]
  molnat[nmol] = fix(s[2])
  for j = 0, molnat[nmol]-1 do begin
    molatid[nmol,j] = s[3+2*j]
    molatmass[nmol,j] = float(s[4+2*j])
  endfor
  nmol = nmol + 1
  skipmol:
endfor

endmolecules:
molid = molid[0:nmol-1]
moldata = moldata[0:nmol-1]
molnat = molnat[0:nmol-1]
molatid = molatid[0:nmol-1, *]
molatmass = molatmass[0:nmol-1, *]
 
free_lun, lun

; calculate thermodynamic properties on the following T grid

;T = 10.^(findgen(101) * 0.06 - 2.)
;T = [0.1, 0.5, 1.0, 3.0, 5.0, 10., 20., 50., 70., findgen(100)*100+100]
;T = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, (findgen(10000)+1.)/10000. * 10000. ]
T = 10.^(findgen(501)/100-1.)
theta = 5040.d0/T
nT = n_elements(T)

; calculate atomic partition functions

Qatom = dblarr(natom, nT)
atom_potion = dblarr(natom)

print, ' '
print, 'Calculating for atoms:'
print, ' '
for i = 0, natom-1 do begin
  print, 'partition function: ',atomid[i]

  read_atom, atomdata[i], J, E, potion
  comp_atom_partf, T, J, E, Qi 
  Qatom[i, *] = Qi
  atom_potion[i] = potion
  
  save, file = atomid[i]+'.sav', Qi, T

endfor  

; below we perform all tasks on each molecule in turn

Qmol = dblarr(nmol, nT)
lgKpmol = dblarr(nmol, nT)
Dmol = dblarr(nmol)

print, ' '
print, 'Calculating for molecules:'
for i = 0, nmol-1 do begin  
  read_mol, moldata[i], molstruct
  Dmol[i] = molstruct.D
  
  print, ' '
  print, molid[i], molstruct.ns, molstruct.D
  print, 'Partition function: ', molid[i]
  comp_mol_partf, T, molstruct, Qi
  Qmol[i, *] = Qi
  
  print, 'Equilibrium constant: ', molid[i]
  i1 = where(atomid eq molatid[i,0], n1) 
  i2 = where(atomid eq molatid[i,1], n2) 
  Qa1 = Qatom[i1,*]
  Qa2 = Qatom[i2,*]
  mu = (molatmass[i,0] * molatmass[i,1]) / (molatmass[i,0] + molatmass[i,1])
  comp_mol_equilkp, T, molstruct, Qa1, Qa2, Qi, mu, lgKp
  lgKpmol[i, *] = lgKp
  
  save, file = molid[i]+'.sav', lgKp, Qi, T
endfor

; save all data in a single file

print, ' '
print, 'saving all data in alldata.sav'
save, file = 'alldata.sav', T, lgKpmol, Qmol, Qatom, molid, atomid, atom_potion, Dmol

; make any comparisons here

; read in Sauval and Tatum 1984 data for comparison

read_sta, aiST, cST, maxE
read_stm, miST, aST, bST, DST
read_ngm, miNG, aNG, bNG, DNG


set_plot, 'ps'
device, file = 'atomic_part.ps', ysize=25
plotthick, 3
!p.charsize=1.2
!p.font=-1
!p.multi=[0,1,2]

print, ' '
print, 'plotting atomic partition functions in atomic_part.ps'

for i = 0, natom-1 do begin 
   plot, T, Qatom[i,*], title = atomid[i], xtitle = 'T [K]', $
                        ytitle = 'Q' 
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, thick=1
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   'NextGen', charsize=0.8 
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, thick=7
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   'Sauval & Tatum (1984)', charsize=0.8			
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 0
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8				

   ;ind = where(aiNG eq atomid[i], nind) 
  ; 
  ; if nind gt 0 then begin
  ;    
  ;    QNG = 10.d0^poly(alog10(theta), cNG[*,ind])
  ;    oplot, T, QNG, linestyle = 2, thick=1 
  ; 
  ; endif   
			
   ind = where(aiST eq atomid[i], nind)

   if nind gt 0 then begin
      
      QST = 10.d0^poly(alog10(theta), cST[*,ind])
      oplot, T, QST, linestyle = 1, thick=7 
      
      ind2 = where(QST lt 200)
      
      yy = (Qatom[i,ind2]-QST[ind2]) / Qatom[i,ind2] * 100.
      yyr = [min(yy)-abs(max(yy)-min(yy))*0.05, $ 
             max(yy)+abs(max(yy)-min(yy))*0.05]
      plot, T[ind2], yy, $
            title = atomid[i], xtitle = 'T [K]', ytitle = '!4D!3Q [%]', $
	    /ys, yr=yyr, thick = 5
      
   endif else begin
      plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4		
   endelse

endfor

device, /close_file
set_plot, 'x'
!p.multi=0


set_plot, 'ps'
device, file = 'molecule_part.ps', ysize=25
plotthick, 3
!p.charsize=2
!p.font=-1
!p.multi=[0,1,3]

print, ' '
print, 'plotting molecular partition functions in molecule_part.ps'

for i = 0, nmol-1 do begin 
   plot, T, Qmol[i,*], title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'Q' 
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, thick=1
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   'NextGen', charsize=0.8
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, thick=7
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   'Sauval & Tatum (1984)', charsize=0.8			
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 0
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8				
	   
   ind = where(miNG eq molid[i], nind)
   
   if nind gt 0 then begin
      QNG = 10.d0^poly(alog10(theta), aNG[*,ind])
      oplot, T, QNG, linestyle = 2, thick=1 
   endif
   		
   ind = where(miST eq molid[i], nind)
   if nind gt 0 then begin
      QST = 10.d0^poly(alog10(theta), aST[*,ind])
      oplot, T, QST, linestyle = 1, thick=7 
      
   endif
   
   plot, T, Qmol[i,*], /ylog, /xlog, title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'Q' 
			
			
			
   ind = where(miNG eq molid[i], nind)
   
   if nind gt 0 then begin
      QNG = 10.d0^poly(alog10(theta), aNG[*,ind])
      oplot, T, QNG, linestyle = 2, thick=1 
   endif			
			
   ind = where(miST eq molid[i], nind)
   
   if nind gt 0 then begin
      QST = 10.d0^poly(alog10(theta), aST[*,ind])
      oplot, T, QST, linestyle = 1, thick=7 
      
      plot, T, (Qmol[i,*]-QST)/Qmol[i,*]*100, /xlog, $
         title = molid[i], xtitle = 'T [K]', $
         ytitle = '!4D!3Q [%]', thick=5, /ys 
   endif else begin
      plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4		 
   endelse
  
   
   
endfor

device, /close_file
set_plot, 'x'
!p.multi=0

set_plot, 'ps'
device, file = 'molecule_equil.ps', ysize=25
plotthick, 3
!p.charsize=2
!p.font=-1
!p.multi=[0,1,3]

print, ' '
print, 'plotting molecular equilibrium functions in molecule_equil.ps'

for i = 0, nmol-1 do begin 
   plot, T, lgKpmol[i,*], title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'log K!dp!n', $
			/xlog, /ys 
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, thick=1
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   'NextGen', charsize=0.8
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, thick=7
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   'Sauval & Tatum (1984)', charsize=0.8			
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 0
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8				
			
			
   ind2 = where(miNG eq molid[i], nind2)
   if nind2 gt 0 then begin      
      lgKpNG = poly(alog10(theta), bNG[*,ind2])-theta*DNG[ind2[0]]
      oplot, T, lgKpNG , psym=4;linestyle=2, thick=1       
   endif
   
   ind = where(miST eq molid[i], nind)
   if nind gt 0 then begin      
      lgKpST = poly(alog10(theta), bST[*,ind])-theta*DST[ind[0]]
      oplot, T, lgKpST , linestyle=1, thick=7       
   endif
   plot, T, lgKpmol[i,*], title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'log K!dp!n', $
			/xlog, xr=[100, 10000], /ys 
   
   if nind2 gt 0 then begin      
      oplot, T, lgKpNG , psym=4;linestyle=2, thick=1       
   endif
   
   if nind gt 0 then begin
      oplot, T, lgKpST , linestyle=1, thick=7        
   
   plot, T, (lgKpmol[i,*]- lgKpST), $
            title = molid[i], xtitle = 'T [K]', $
                        ytitle = '!4D!3 log K!dp!n', $
			/xlog, xr=[100, 10000], thick = 5
   endif else begin
   
      plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4		
   endelse
      		 
   
endfor

device, /close_file
set_plot, 'x'
plotthick, 1
!p.multi=0

print, ' '
print, 'normal end'
print, ' '
stop
end


