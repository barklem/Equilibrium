pro collate

; this routine organises all the data for publication
; 
; Paul Barklem, February 2012

; first 

outdir1 = './Final_Collated_Results/'
outdir2 = './Final_Collated_Results/Mols_ascii/'
indir1 = './PartFuncs_Results/'
indir2 = './Final_Mols_ascii/'


; collate dissociation energies
read_diss, molname1, De, hh=hh, luo=luo, g2=g2, other=other, components=molcomp
name = molname1
restore, indir1 + 'allpartf.idl'  ; gives : T, Qmol, Qmol4, QmolSt, QmolHH, molid, Dmol, EZEROmol


; clean out cases with no data
nmol = n_elements(Dmol)
tqmol = dblarr(nmol)
for i = 0, nmol - 1 do tqmol[i] = total(Qmol[i, *])
ind = where(tqmol gt 0.d0, nind, complement = ind2)
indstore = ind
nindstore = nind
Qmol = Qmol[ind,*]
Qmol4 = Qmol4[ind,*]
QmolSt = QmolSt[ind,*]
QmolHH = QmolHH[ind,*]
molid = molid[ind]
EZEROmol = EZEROmol[ind]
molid_cutout = molid[ind2]
nmol = nind

; don't assume molecule lists are identical, but they should be
name_new = strarr(nmol)
molcomp_new = strarr(nmol,2)
hh_new = dblarr(nmol,2)
luo_new = dblarr(nmol,2)
g2_new = dblarr(nmol,2)
other_new = dblarr(nmol,2)
De_new = dblarr(nmol,2)

hh_s = strarr(nmol,2)+'.'
luo_s = strarr(nmol,2)+'.'
g2_s = strarr(nmol,2)+'.'
other_s = strarr(nmol,2)+'.'
De_s = strarr(nmol,2)+'.'

for i = 0, nmol-1 do begin
ind3 = where(molname1 eq molid[i], nind3)
  if nind3 le 0 then stop

name_new[i] = name[ind3]
molcomp_new[i,*] = molcomp[ind3,*]
hh_new[i,*] = hh[ind3, *]
luo_new[i,*] = luo[ind3,*]
g2_new[i,*] = g2[ind3,*]
other_new[i,*] = other[ind3,*]
De_new[i,*] = De[ind3,*]

for j = 0, 1 do begin
if hh_new[i,j] gt 0.d0 then hh_s[i,j] = string(hh[ind3, j], '(f10.6)')
if luo_new[i,j] gt 0.d0 then luo_s[i,j] = string(luo[ind3, j], '(f10.6)')
if g2_new[i,j] gt 0.d0 then g2_s[i,j] = string(g2[ind3, j], '(f10.6)')
if other_new[i,j] gt 0.d0 then other_s[i,j] = string(other[ind3, j], '(f10.6)')
if De_new[i,j] gt 0.d0 then De_s[i,j] = string(De[ind3, j], '(f10.6)')
endfor

endfor

name = name_new
molcomp = molcomp_new
hh = hh_new
luo = luo_new
g2 = g2_new
other = other_new
De = De_new





spawn, 'rm '+outdir2+'*.txt'

openw, 1, outdir1 + 'diss_table.txt'
printf, 1, 'Dissociation energies', nmol, ' molecules
printf, 1, 'Molecule', 'HH', 'Luo', 'G2', 'Adopted', format = '(a15, 5a22)
printf, 1, ' ', 'Do', 'sigma(Do)', 'Do', 'sigma(Do)', 'Do', 'sigma(Do)', 'Do', 'sigma(Do)', format = '(a15, 10a11)
printf, 1, ' ', '[eV]', '[eV]', '[eV]', '[eV]', '[eV]', '[eV]', '[eV]', '[eV]', format = '(a15, 10a11)

for i = 0, nmol - 1 do begin
      printf, 1, molid[i], molcomp[i,0], molcomp[i,1], hh_s[i,0], hh_s[i,1], luo_s[i,0], luo_s[i,1], g2_s[i,0], g2_s[i,1], $ 
      ; other_s[i,0], other_s[i,1], 
      De_s[i,0], De_s[i,1], format = '(3a5,5(2x,2a10))'
      spawn, 'cp '+indir2+molid[i]+'.txt '+outdir2+molid[i]+'.txt'
endfor

close, 1


openw, 1, outdir1 + 'partf_table.txt'
printf, 1, 'Partition functions Q'
printf, 1, nind
printf, 1, T, format = '("T [K]", 600(2x,e12.5))'
printf, 1, ' '
for i = 0, nind-1 do printf, 1, molid[i], Qmol[i, *], format = '(a5, 600(2x,e12.5))'
close, 1

openw, 1, outdir1 + 'ezero_table.txt'
printf, 1, 'Energy of first state wrt bottom of potential well: Ezero'
printf, 1, nind 
printf, 1, ' '
for i = 0, nind-1 do printf, 1, molid[i], EZEROmol[i], format = '(a5, 2x,e12.5)'
close, 1

; new plots only including these



set_plot, 'ps'
device, file = outdir1+'diss.eps', xsize=30, ysize=18
plotthick, 3
!p.charsize = 2.0

!p.multi = [0,3,2]

psymcircle
ind = where(luo(*,0) gt 0.01 and hh(*,0) gt 0.01, nind)
plot, luo(ind,0), hh(ind,0), psym = 8, ytitle = 'D!d0!n!u0!n(HH) [eV]', xtitle = 'D!d0!n!u0!n(Luo) [eV]'
oplot, [0,20], [0,20], linestyle = 1
for i = 0, nind-1 do begin
   if luo(ind(i),1) ne 0.d0 then begin
       oplot, [luo(ind(i),0)-luo(ind(i),1), luo(ind(i),0)+luo(ind(i),1)], [hh(ind(i),0), hh(ind(i),0)], thick=2
   endif
endfor
for i = 0, nind-1 do begin
   if hh(ind(i),1) ne 0.d0 then begin
       oplot, [luo(ind(i),0), luo(ind(i),0)], [hh(ind(i),0)-hh(ind(i),1), hh(ind(i),0)+hh(ind(i),1)], thick=2
   endif
endfor

ind = where(luo(*,0) gt 0.01 and g2(*,0) gt 0.01, nind)
plot, luo(ind,0), g2(ind,0), psym = 8, ytitle = 'D!d0!n!u0!n(G2) [eV]', xtitle = 'D!d0!n!u0!n(Luo) [eV]' 
oplot, [0,20], [0,20], linestyle = 1
for i = 0, nind-1 do begin
   if luo(ind(i),1) ne 0.d0 then begin
       oplot, [luo(ind(i),0)-luo(ind(i),1), luo(ind(i),0)+luo(ind(i),1)], [g2(ind(i),0), g2(ind(i),0)], thick=2
   endif
endfor

ind = where(g2(*,0) gt 0.01 and hh(*,0) gt 0.01, nind)
plot, g2(ind,0), hh(ind,0), psym = 8, xtitle = 'D!d0!n!u0!n(G2) [eV]', ytitle = 'D!d0!n!u0!n(HH) [eV]'
oplot, [0,20], [0,20], linestyle = 1
for i = 0, nind-1 do begin
   if hh(ind(i),1) ne 0.d0 then begin
       oplot, [g2(ind(i),0), g2(ind(i),0)], [hh(ind(i),0)-hh(ind(i),1), hh(ind(i),0)+hh(ind(i),1)], thick=2
   endif
endfor

openw, 1, outdir1 + 'diss_luo_HH.txt'
printf, 1, 'Luo vs HH'
printf, 1, '        Luo    HH     diff'
diff = hh(*,0) - luo(*,0) 
diff1 = abs(diff)
ind2 = where(hh(*,0) ne 0.d0 and luo(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
printf, 1, name(ind2(ind(i))), luo(ind2(ind(i)),0), hh(ind2(ind(i)),0), diff(ind2(ind(i))), diff(ind2(ind(i)))/luo(ind2(ind(i)),0), format = '(a5, 4f7.2)'
endfor
close, 1

plot, luo(ind2,0), diff(ind2), psym = 8, /ys, yr=[-2, 2], xtitle = 'D!d0!n!u0!n(Luo) [eV]', ytitle = '!4D!XD!d0!n!u0!n [eV]'

for i = 0, nind2-1 do begin
   oplot, [luo(ind2[i],0), luo(ind2[i],0)], [diff(ind2[i]) - luo(ind2[i],1) - hh(ind2[i],1), diff(ind2[i]) + luo(ind2[i],1) + hh(ind2[i],1)], thick=2
endfor


openw, 1, outdir1 + 'diss_luo_g2.txt'
printf, 1, 'Luo vs G2'
printf, 1, '        Luo    G2     diff'
diff =  g2(*,0) - luo(*,0) 
diff1 = abs(diff)
ind2 = where(g2(*,0) ne 0.d0 and luo(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
printf, 1, name(ind2(ind(i))), luo(ind2(ind(i)),0), g2(ind2(ind(i)),0), diff(ind2(ind(i))),  diff(ind2(ind(i)))/luo(ind2(ind(i)),0), format = '(a5, 4f7.2)'
endfor
close, 1

plot, luo(ind2,0), diff(ind2), psym = 8, /ys, yr=[-2, 2], xtitle = 'D!d0!n!u0!n(Luo) [eV]', ytitle = '!4D!XD!d0!n!u0!n [eV]'

for i = 0, nind2-1 do begin
   oplot, [luo(ind2[i],0), luo(ind2[i],0)], [diff(ind2[i]) - luo(ind2[i],1) - g2(ind2[i],1), diff(ind2[i]) + luo(ind2[i],1) + g2(ind2[i],1)], thick=2
endfor


openw, 1, outdir1 + 'diss_HH_g2.txt'
printf, 1, 'HH vs G2'
printf, 1, '        HH     G2     diff'
diff = hh(*,0) - g2(*,0)
diff1 = abs(diff)
ind2 = where(g2(*,0) ne 0.d0 and hh(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
printf, 1, name(ind2(ind(i))), hh(ind2(ind(i)),0), g2(ind2(ind(i)),0), diff(ind2(ind(i))),  diff(ind2(ind(i)))/luo(ind2(ind(i)),0), format = '(a5, 4f7.2)'
endfor
close, 1

plot, g2(ind2,0), diff(ind2), psym = 8, /ys, yr=[-2, 2], xtitle = 'D!d0!n!u0!n(G2) [eV]', ytitle = '!4D!XD!d0!n!u0!n [eV]'

for i = 0, nind2-1 do begin
   oplot, [g2(ind2[i],0), g2(ind2[i],0)], [diff(ind2[i]) - g2(ind2[i],1) - hh(ind2[i],1), diff(ind2[i]) + g2(ind2[i],1) + hh(ind2[i],1)], thick=2
endfor

device, /close_file
set_plot, 'x'



set_plot, 'ps'
device, file = outdir1+'diss2.eps', /portrait, ysize=30
plotthick, 3
!p.charsize = 1.3

!p.multi = [0,1,2]

psymcircle
ind = where(De(*,0) gt 0.01 and hh(*,0) gt 0.01, nind)
plot, De(ind,0), hh(ind,0), psym = 8, ytitle = 'D!d0!n!u0!n(HH) [eV]', xtitle = 'D!d0!n!u0!n(Adopted value) [eV]'
oplot, [0,20], [0,20], linestyle = 1

for i = 0, nind-1 do begin
   if hh(ind(i),1) ne 0.d0 then begin
       oplot, [De(ind(i),0), De(ind(i),0)], [hh(ind(i),0)-hh(ind(i),1), hh(ind(i),0)+hh(ind(i),1)], thick=2
   endif
endfor

for i = 0, nind-1 do begin
   if De(ind(i),1) ne 0.d0 then begin
       oplot, [De(ind(i),0)-De(ind(i),1), De(ind(i),0)+De(ind(i),1)], [hh(ind(i),0), hh(ind(i),0)], thick=2
   endif
endfor

openw, 1, outdir1 + 'diss_adopted_HH.txt'
printf, 1, 'Dissociation energies - Differences between HH and adopted values'
printf, 1, 'Ordered by absolute difference, largest to smallest'
printf, 1, '        Luo    HH     diff   rel-diff'
diff = hh(*,0) - De(*,0) 
diff1 = abs(diff)
ind2 = where(hh(*,0) ne 0.d0 and De(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
printf, 1, name(ind2(ind(i))), De(ind2(ind(i)),0), hh(ind2(ind(i)),0), diff(ind2(ind(i))), diff(ind2(ind(i)))/De(ind2(ind(i))), format = '(a5, 4f7.2)'
endfor
close, 1

ind3 = ind2[ind[0:12]]

plot, De(ind2,0), diff(ind2), psym = 8, /ys, yr=[-2, 2], xtitle = 'D!d0!n!u0!n(Adopted value) [eV]', ytitle = '!4D!XD!d0!n!u0!n [eV]', symsize = 1.0

for i = 0, nind2-1 do begin
   oplot, [De(ind2[i],0), De(ind2[i],0)], [diff(ind2[i]) - De(ind2[i],1), diff(ind2[i]) + De(ind2[i],1) ], thick=2
endfor
xyouts, De[ind3,0]+0.1, diff[ind3], name[ind3], charsize = 1.0

device, /close_file
set_plot, 'x'



read_stm, miST, aST, bST, DST

nT = n_elements(T)


Tc = [1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000.]
theta = 5040.d0/Tc

inds = [where(molid eq 'H2'), where(molid eq 'CO')]
for j = 0, n_elements(inds)-1 do begin
   i = inds[j]


set_plot, 'ps'
device, file = outdir1 + molid[i]+'_partf.eps', /color;, ysize=25, yoffset = 2
colors
plotthick, 3
!p.charsize=2
!p.font=-1
!p.multi=[0,1,3]

	   
   		
   ind = where(miST eq molid[i], nind)
   if nind gt 0 then begin
      QST = 10.d0^poly(alog10(theta), aST[*,ind])
      
   endif
   
   ; Irwin 1987
   if molid[i] eq 'H2' then begin
      AI87 = [1.67298118410d4, -1.49945289142d4, 5.74838863349d3, -1.22210505066d3, $
             1.55637569965d2, -1.18744926193d1, 5.02617615447d-1,-9.10563051348d-3]
      QI87 = exp(poly(alog(Tc), AI87))
   endif
   if molid[i] eq 'CO' then begin
      AI87 = [5.05610415417d4, -5.19600025580d4, 2.33277267148d4, -5.97599449706d3, $
              9.55509531681d2, -9.76517012179d1, 6.22988547018d0, -2.26856284960d-1, $
              3.61025385248d-3]
      QI87 = exp(poly(alog(Tc), AI87))
   endif
   
   QST = T*0.
   ind2 = where(miST eq molid[i], nind2)
   if nind2 gt 0 then begin
      QST = 10.d0^poly(alog10(theta), aST[*,ind2])
      ymax = max([transpose(Qmol[i,*]), QST, 1.])
      ymin = min([transpose(Qmol[i,*]), QST, 0.1])
   endif else begin
      ymax = max([transpose(Qmol[i,*]), 1.])
      ymin = min([transpose(Qmol[i,*]), 0.1])
   endelse
   ymin = max([ymin, 1e-20])
   ;diff = alog10(ymax) - alog10(ymin)
   ;ymax = 10.^(alog10(ymax) + 0.05 * diff)
   ;ymin = 10.^(alog10(ymin) - 0.05 * diff)
   ;print, ymin, ymax, diff


Tg = 0.
  if molid[i] eq 'CO' then begin
      GGL = [0.27758,   0.36290, -0.74669e-5, 0.14896e-7]
      GGM = [0.90723e1, 0.33263,  0.11806e-4, 0.27035e-7]
      GGH = [0.63418e2, 0.20760,  0.10895e-3, 0.19844e-8]
      QL = poly(T, GGL)
      QM = poly(T, GGM)
      QH = poly(T, GGH)
      igl = where((T ge 70) and (T le 500))
      igm = where((T gt 500) and (T le 1500))
      igh = where((T gt 1500) and (T le 3005))
      Tg = [T[igl], T[igm], T[igh]]
      QGG = [QL[igl], QM[igm], QH[igh]]
   endif 

    plot, T, Qmol[i,*], /ylog, /xlog, xtitle = 'T [K]', $
                        ytitle = 'Q', yr = [ymin, ymax], linestyle = 1, thick = 3
                    			
			
   
   if molid[i] eq 'H2' or molid[i] eq 'CO' then oplot, Tc, QI87, linestyle = 5, color =2, thick = 3, psym = -1
			
   
   
   if nind2 gt 0 then begin      
      oplot, Tc, QST, linestyle = 2, color = 4, thick = 3, psym = -6
      
     ; if total(Qmol[i,*]) gt 0.d0 then $
     ; plot, T, (Qmol[i,*]-QST)/Qmol[i,*]*100, /xlog, $
     ;    title = molid[i], xtitle = 'T [K]', $
     ;    ytitle = '!4D!3Q [%] (This work - ST)', thick=3, /ys 
   endif ;else begin
     ; plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4		 
   ;endelse

   
   
   plot, T, Qmol[i,*], xtitle = 'T [K]', $
                        ytitle = 'Q', /nodata
oplot, T, Qmol[i, *], linestyle = 1, thick =3
     
 oplot, [!x.crange[0] + 0.1*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.2*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.70*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.70*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 5, color = 2, thick = 3, psym = -1
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.70*(!y.crange[1]-!y.crange[0]), $
	   'Irwin (1987)', charsize=0.8 
   oplot, [!x.crange[0] + 0.1*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.2*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.60*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.60*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, color = 4, thick = 3, psym = -6
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.60*(!y.crange[1]-!y.crange[0]), $
	   'Sauval & Tatum (1984)', charsize=0.8			
   oplot, [!x.crange[0] + 0.1*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.2*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.50*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.50*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, thick = 3
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.5*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8	
	      
	   
   		
   ind = where(miST eq molid[i], nind)
   if nind gt 0 then begin
      QST = 10.d0^poly(alog10(theta), aST[*,ind])
      oplot, Tc, QST, linestyle = 2, color = 4, thick = 3, psym =-6
      
   endif
   
   ; Irwin 1987
   if molid[i] eq 'H2' then begin
      oplot, Tc, QI87, linestyle = 5, color = 2, thick = 3, psym = -1
   endif
   if molid[i] eq 'CO' then begin
      oplot, Tc, QI87, linestyle = 5, color = 2, thick = 3, psym = -1
   endif
   
   

 plot, T, Qmol[i,*], xtitle = 'T [K]', $
                        ytitle = 'Q' , /nodata
   oplot,  T, Qmol[i,*], linestyle = 1, thick = 5
   oplot, T, QmolHH[i,*], linestyle = 5, color = 2, thick = 3
   oplot, T, Qmol4[i,*], linestyle = 3, color = 3, thick = 3
   oplot, T, QmolST[i,*], linestyle = 2, color = 4, thick = 3
   
   QmHH = interpol(QmolHH[i,*], T, Tc)
   Qm4 = interpol(Qmol4[i,*], T, Tc)
   QmSt = interpol(QmolSt[i,*], T, Tc)
   oplot, Tc, QmHH,  color = 2, thick = 3, psym = 1
   oplot, Tc, Qm4,  color = 3, thick = 3, psym = 4
   oplot, Tc, QmST,  color = 4, thick = 3, psym = 6
   
   oplot, [!x.crange[0] + 0.10*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.20*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.70*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.70*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 5, color = 2, thick = 3, psym = -1
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.70*(!y.crange[1]-!y.crange[0]), $
	   'HH', charsize=0.8
   oplot, [!x.crange[0] + 0.10*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.20*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.60*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.60*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, color = 4 , thick = 3, psym = -6
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.60*(!y.crange[1]-!y.crange[0]), $
	   'HH + no high order (=ST)', charsize=0.8
   oplot, [!x.crange[0] + 0.10*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.20*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.50*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.50*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 3, color = 3, thick = 3, psym=-4
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.50*(!y.crange[1]-!y.crange[0]), $
	   'States < 40000 /cm', charsize=0.8			
   oplot, [!x.crange[0] + 0.10*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.20*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.40*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.40*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, thick = 5
   xyouts, !x.crange[0] + 0.25*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.40*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8

   
   device, /close
   set_plot, 'x'
   !p.multi=0

; DO ALTERNATE PLOT IN REMO'S STYLE


cd, current=current_dir

; plot parameters ps files
encapsulated  = 1
landscape     = 0
font          = -1
default_thick = 4
charsize      = 1.2
charthick     = default_thick

; plot parameters, other
thick         = 6
thin          = 3
thick_nist    = 9 
thick_irwin   = 7 
thick_moog    = 7
thick_extra    = 7
xtitle        = 'Temperature [K]'
ytitle_top    = 'Q (partition function)'
ytitle_bot    = '$\Delta$Q / Q'
color_nist    = 'blu4'
color_irwin   = 'tomato'
color_moog    = 'grn4'
color_extra   = 'plum'
linest_nist   = '0'
linest_irwin  = '5'
linest_moog   = '2'
linest_extra  = '3'

psym_moog     = 'filled circle'
syms_moog    = 1.2

; plot: x and y ranges
xmin          =  1.0e-2
xmax          =  1.0e+4
ymin_bot      = -0.23
ymax_bot      = +0.08
      

Qmol_xirwin = exp(interpol(alog(Qmol[i,*]),alog(T),alog(Tc),/spline))
Qmol_xST = exp(interpol(alog(Qmol[i,*]),alog(T),alog(Tc),/spline))
Qmol_xG = exp(interpol(alog(Qmol[i,*]),alog(T),alog(Tg),/spline))

; plot
; open PS file
      cgPS_open, outdir1 + molid[i]+'_partf2.eps', encapsulated=encapsulated, landscape=landscape, $
                 font=font, charsize=charsize, default_thick=default_thick, $
                 /nomatch
; plot top panel
; y range, min and max values
;;      ymin_top = min([reform(qatom_nist[j_nist,*]),reform(qatom_irwin),reform(qmoog_tnist)], max=ymax_top)
      ymin_top = min([reform(Qmol[i,*])], max=ymax_top)
; extend range slightly
      ymin_top = ymin_top*0.5
      ymax_top = ymax_top*2.0
; custom tickmarks, logarithmic y axis
      ticks  = LogLevels([ymin_top,ymax_top],/fine)
      nticks = N_Elements(ticks)
      cgplot,    T, Qmol[i,*], linest=linest_nist, thick=thick_nist, col=color_nist, $
                 XTickformat='(A1)', ytitle=ytitle_top,  title=molid[i], $
                 /xlog, XRange=[xmin,xmax], $
                 /ylog, YRange=[ymin_top,ymax_top], YStyle=2, YTicks=nticks-1, YTickV=ticks, $
                 position=[0.17,0.32,0.90,0.90]
      if molid[i] eq 'H2' or molid[i] eq 'CO' then cgplot, /ov, Tc, QI87, linest=linest_irwin, thick=thick_irwin, col=color_irwin
      if nind2 gt 0 then cgplot, /ov, Tc, QST, linest=linest_moog, thick=thick_moog, col=color_moog
      if molid[i] eq 'CO' then cgplot, /ov, Tg, Qgg, linest=linest_extra, thick=thick_extra, col=color_extra
      
      ;ww_irwin = where(temp_nist ge min_temp_irwin)
      ;cww_irwin= where(temp_nist le min_temp_irwin)
      ;cgplot,/ov,temp_nist[ ww_irwin], qatom_irwin[ ww_irwin], linest=linest_nist, thick=thick_irwin, col=color_irwin
      ;cgplot,/ov,temp_nist[cww_irwin], qatom_irwin[cww_irwin], linest=linest_irwin, thick=thick_irwin, col=color_irwin
      ;cgplot,/ov,tmoog, qmoog, psym=cgsymcat(psym_moog),col=color_moog,syms=syms_moog
      ;ww_moog = where(temp_nist ge min_tmoog)
      ;cww_moog= where(temp_nist le min_tmoog)
      ;cgplot,/ov,temp_nist, qmoog_tnist, linest=linest_moog, thick=thick_moog, col=color_moog

; legend

      if molid[i] eq 'CO' then begin
      xloc   = 10.^(!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.05)
      yloc   = 10.^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.9)
      cgLegend, Title=['This work', 'Sauval & Tatum (1984)', 'Irwin (1987)', 'Gamache et al. (2000)'], $
                Lines=[linest_nist,linest_moog,linest_irwin,linest_extra], $
                ;PSym=[0,0,cgsymcat(psym_moog)], syms=syms_moog, $ ;syms=[1,1,syms_moog], $
                thick=thick_nist, $
                Color=[color_nist,color_moog,color_irwin,color_extra], $
                charsize=charsize, charthick=charthick, vspace=1.8, length=0.125, $
                Location=[xloc,yloc], /Data, /background

      endif else begin
      xloc   = 10.^(!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.05)
      yloc   = 10.^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.9)
      cgLegend, Title=['This work', 'Sauval & Tatum (1984)', 'Irwin (1987)'], $
                Lines=[linest_nist,linest_moog,linest_irwin], $
                ;PSym=[0,0,cgsymcat(psym_moog)], syms=syms_moog, $ ;syms=[1,1,syms_moog], $
                thick=thick_nist, $
                Color=[color_nist,color_moog,color_irwin], $
                charsize=charsize, charthick=charthick, vspace=1.8, length=0.125, $
                Location=[xloc,yloc], /Data, /background
      endelse

; plot bottom panel
      cgplot,/noerase, [xmin,xmax], [0,0], linest=linest_nist, thick=thick_nist, col=color_nist, $
             xtitle=xtitle, ytitle=ytitle_bot, $
             /xlog, xrange=[xmin,xmax], xstyle=8, $
             yrange=[ymin_bot,ymax_bot], ystyle=1, $
             position=[0.17,0.15,0.90,0.32], xticklen=0.1
      cgplot,/ov,Tc,(QI87/Qmol_xirwin-1.), linest=linest_irwin, thick=thick_irwin, col=color_irwin
      cgplot,/ov,Tc,(QST/Qmol_xST-1.), linest=linest_moog, thick=thick_moog, col=color_moog
      if molid[i] eq 'CO' then cgplot,/ov,Tg,(Qgg/Qmol_xG-1.), linest=linest_extra, thick=thick_extra, col=color_extra


; close PS file
cgPS_close
  
xmin = 1000.

; plot
; open PS file
      cgPS_open, outdir1 + molid[i]+'_partf3.eps', encapsulated=encapsulated, landscape=landscape, $
                 font=font, charsize=charsize, default_thick=default_thick, $
                 /nomatch
; plot top panel
; y range, min and max values
;;      ymin_top = min([reform(qatom_nist[j_nist,*]),reform(qatom_irwin),reform(qmoog_tnist)], max=ymax_top)
      ymin_top = min([reform(Qmol[i,*])], max=ymax_top)
; extend range slightly
      ymin_top = ymin_top*0.5
      ymax_top = ymax_top*1.1
; custom tickmarks, logarithmic y axis
      ticks  = LogLevels([ymin_top,ymax_top],/fine)
      nticks = N_Elements(ticks)
      cgplot,    T, Qmol[i,*], linest=linest_nist, thick=thick_nist, col=color_nist, $
                 XTickformat='(A1)', ytitle=ytitle_top,  title=molid[i], $
                 XRange=[xmin,xmax], $
                 YRange=[0.,ymax_top], $;YStyle=2, YTicks=nticks-1, YTickV=ticks, $
                 position=[0.17,0.32,0.90,0.90]
      if molid[i] eq 'H2' or molid[i] eq 'CO' then cgplot, /ov, Tc, QI87, linest=linest_irwin, thick=thick_irwin, col=color_irwin
      if nind2 gt 0 then cgplot, /ov, Tc, QST, linest=linest_moog, thick=thick_moog, col=color_moog
      if molid[i] eq 'CO' then cgplot, /ov, Tg, Qgg, linest=linest_extra, thick=thick_extra, col=color_extra
      
      ;ww_irwin = where(temp_nist ge min_temp_irwin)
      ;cww_irwin= where(temp_nist le min_temp_irwin)
      ;cgplot,/ov,temp_nist[ ww_irwin], qatom_irwin[ ww_irwin], linest=linest_nist, thick=thick_irwin, col=color_irwin
      ;cgplot,/ov,temp_nist[cww_irwin], qatom_irwin[cww_irwin], linest=linest_irwin, thick=thick_irwin, col=color_irwin
      ;cgplot,/ov,tmoog, qmoog, psym=cgsymcat(psym_moog),col=color_moog,syms=syms_moog
      ;ww_moog = where(temp_nist ge min_tmoog)
      ;cww_moog= where(temp_nist le min_tmoog)
      ;cgplot,/ov,temp_nist, qmoog_tnist, linest=linest_moog, thick=thick_moog, col=color_moog

; legend


      if molid[i] eq 'CO' then begin
      xloc   = 10.^(!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.05)
      yloc   = 10.^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.9)
      cgLegend, Title=['This work', 'Sauval & Tatum (1984)', 'Irwin (1987)', 'Gamache et al. (2000)'], $
                Lines=[linest_nist,linest_moog,linest_irwin,linest_extra], $
                ;PSym=[0,0,cgsymcat(psym_moog)], syms=syms_moog, $ ;syms=[1,1,syms_moog], $
                thick=thick_nist, $
                Color=[color_nist,color_moog,color_irwin,color_extra], $
                charsize=charsize, charthick=charthick, vspace=1.8, length=0.125, $
                Location=[xloc,yloc], /Data, /background

      endif else begin
      xloc   = 10.^(!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.05)
      yloc   = 10.^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.9)
      cgLegend, Title=['This work', 'Sauval & Tatum (1984)', 'Irwin (1987)'], $
                Lines=[linest_nist,linest_moog,linest_irwin], $
                ;PSym=[0,0,cgsymcat(psym_moog)], syms=syms_moog, $ ;syms=[1,1,syms_moog], $
                thick=thick_nist, $
                Color=[color_nist,color_moog,color_irwin], $
                charsize=charsize, charthick=charthick, vspace=1.8, length=0.125, $
                Location=[xloc,yloc], /Data, /background
      endelse

; plot bottom panel
      cgplot,/noerase, [xmin,xmax], [0,0], linest=linest_nist, thick=thick_nist, col=color_nist, $
             xtitle=xtitle, ytitle=ytitle_bot, $
             /xlog, xrange=[xmin,xmax], xstyle=8, $
             yrange=[ymin_bot,ymax_bot], ystyle=1, $
             position=[0.17,0.15,0.90,0.32], xticklen=0.1
      cgplot,/ov,Tc,(QI87/Qmol_xirwin-1.), linest=linest_irwin, thick=thick_irwin, col=color_irwin
      cgplot,/ov,Tc,(QST/Qmol_xST-1.), linest=linest_moog, thick=thick_moog, col=color_moog
      if molid[i] eq 'CO' then cgplot,/ov,Tg,(Qgg/Qmol_xG-1.), linest=linest_extra, thick=thick_extra, col=color_extra


; close PS file
cgPS_close
  
xmin = 1000.

; plot
; open PS file
      cgPS_open, outdir1 + molid[i]+'_partf4.eps', encapsulated=encapsulated, landscape=landscape, $
                 font=font, charsize=charsize, default_thick=default_thick, $
                 /nomatch
; plot top panel
; y range, min and max values
;;      ymin_top = min([reform(qatom_nist[j_nist,*]),reform(qatom_irwin),reform(qmoog_tnist)], max=ymax_top)
      ymin_top = min([reform(Qmol[i,*])], max=ymax_top)
; extend range slightly
      ymin_top = ymin_top*0.5
      ymax_top = ymax_top*1.1
; custom tickmarks, logarithmic y axis
      ticks  = LogLevels([ymin_top,ymax_top],/fine)
      nticks = N_Elements(ticks)
      cgplot,    T, Qmol[i,*], linest=linest_nist, thick=thick_nist, col=color_nist, $
                 XTickformat='(A1)', ytitle=ytitle_top,  title=molid[i], $
                 XRange=[xmin,xmax], $
                 YRange=[0.,ymax_top], $;YStyle=2, YTicks=nticks-1, YTickV=ticks, $
                 position=[0.17,0.32,0.90,0.90]
      cgplot, /ov, T, QmolHH[i,*], linest=linest_irwin, thick=thick_irwin, col=color_irwin
      cgplot, /ov, T, Qmol4[i,*], linest=linest_extra, thick=thick_moog, col=color_extra
      cgplot, /ov, T, QmolST[i,*], linest=linest_moog, thick=thick_moog, col=color_moog
      
      ;ww_irwin = where(temp_nist ge min_temp_irwin)
      ;cww_irwin= where(temp_nist le min_temp_irwin)
      ;cgplot,/ov,temp_nist[ ww_irwin], qatom_irwin[ ww_irwin], linest=linest_nist, thick=thick_irwin, col=color_irwin
      ;cgplot,/ov,temp_nist[cww_irwin], qatom_irwin[cww_irwin], linest=linest_irwin, thick=thick_irwin, col=color_irwin
      ;cgplot,/ov,tmoog, qmoog, psym=cgsymcat(psym_moog),col=color_moog,syms=syms_moog
      ;ww_moog = where(temp_nist ge min_tmoog)
      ;cww_moog= where(temp_nist le min_tmoog)
      ;cgplot,/ov,temp_nist, qmoog_tnist, linest=linest_moog, thick=thick_moog, col=color_moog

; legend
      xloc   = (!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.05)
      yloc   = (!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.9)
      cgLegend, Title=['This work', 'HH', 'States < 40000 cm!e-1!n', 'HH + no high order (=ST)' ], $
                Lines=[linest_nist,linest_irwin,linest_extra,linest_moog], $
                ;PSym=[0,0,cgsymcat(psym_moog)], syms=syms_moog, $ ;syms=[1,1,syms_moog], $
                thick=thick_nist, $
                Color=[color_nist,color_irwin,color_extra,color_moog], $
                charsize=charsize, charthick=charthick, vspace=1.8, length=0.125, $
                Location=[xloc,yloc], /Data, /background

; plot bottom panel
      cgplot,/noerase, [xmin,xmax], [0,0], linest=linest_nist, thick=thick_nist, col=color_nist, $
             xtitle=xtitle, ytitle=ytitle_bot, $
             /xlog, xrange=[xmin,xmax], xstyle=8, $
             yrange=[ymin_bot,ymax_bot], ystyle=1, $
             position=[0.17,0.15,0.90,0.32], xticklen=0.1
      cgplot,/ov,T,(QmolHH[i,*]/Qmol[i,*]-1.), linest=linest_irwin, thick=thick_irwin, col=color_irwin
      cgplot,/ov,T,(Qmol4[i,*]/Qmol[i,*]-1.), linest=linest_extra, thick=thick_moog, col=color_extra
      cgplot,/ov,T,(QmolST[i,*]/Qmol[i,*]-1.), linest=linest_moog, thick=thick_moog, col=color_moog


; close PS file
cgPS_close


  
endfor

; convert EPS files to PDF
cd, outdir1
spawn,'for i in *'+'.eps'+'; do epstopdf $i ; done'
cd, current_dir


; now output equilibrium constants

restore, indir1 + 'allequil.idl'  ; gives : T, lgKpmol, molid, Dmol, molecule_components

molid = molid[indstore]
lgKpmol = lgKpmol[indstore,*]
Dmol = Dmol[indstore]
molecule_components = molecule_components[indstore,*]

openw, 1, outdir1 + 'equil_table.txt'
printf, 1, 'Equilibrium constants  log10(pK)'
printf, 1, nindstore
printf, 1, T, format = '("T [K]", 600(2x,e12.5))'
printf, 1, ' '
for i = 0, nindstore-1 do printf, 1, molid[i], lgKpmol[i, *], format = '(a5, 600(2x,e12.5))'
close, 1





; now output atomic partition functions

restore, indir1 + 'allatom.idl'  ; gives : T, Qatom, atomid, atom_potion

na = n_elements(atomid)

openw, 1, outdir1 + 'atompartf_table.txt'
printf, 1, 'Partition functions Q'
printf, 1, na
printf, 1, T, format = '("  T [K]", 600(2x,e12.5))'
printf, 1, ' '
for i = 0, na-1 do printf, 1, atomid[i], Qatom[i, *], format = '(a7, 600(2x,e12.5))'
close, 1



save, file = outdir1 + 'allatom_collated.idl', T, Qatom, atomid, atom_potion
save, file = outdir1 + 'allpartf_collated.idl', T, Qmol, Qmol4, QmolSt, QmolHH, molid, Dmol, EZEROmol
save, file = outdir1 + 'allequil_collated.idl', T, lgKpmol, molid, Dmol, molecule_components




stop





end

