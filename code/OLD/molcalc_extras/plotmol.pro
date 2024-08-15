pro plotmol

; file locations 
moldir = './Final_Mols_idl/'  ; should have ending slash
outdir1 = './PartFuncs_Results/'
restore, outdir1 + 'allpartf.idl'

nmol = n_elements(molid)
theta = 5040.d0/T
nT = n_elements(T)

; make any comparisons here

; read in Sauval and Tatum 1984 data for comparison

;read_sta, aiST, cST, maxE
read_stm, miST, aST, bST, DST
read_ngm, miNG, aNG, bNG, DNG

goto, skipatom

set_plot, 'ps'
device, file = 'atomic_part.ps', ysize=25, yoffset=2
plotthick, 3
!p.charsize=1.2
!p.font=0
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

skipatom:

set_plot, 'ps'
device, file = outdir1 + 'allpartf.eps', /landscape, /color;, ysize=25, yoffset = 2
colors
plotthick, 3
!p.charsize=1
!p.font=-1
!p.multi=[0,2,2]

print, ' '
print, 'plotting molecular partition functions in allpartf.ps'

for i = 0, nmol-1 do begin 
   print, molid[i]
   plot, T, Qmol[i,*], title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'Q' 
   oplot,  T, Qmol[i,*], thick = 3
   oplot, T, QmolHH[i,*], linestyle = 5, color = 4
   oplot, T, Qmol4[i,*], linestyle = 2, color = 2
   oplot, T, QmolST[i,*], linestyle = 1, color = 3
   
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 5, color = 4
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   'HH', charsize=0.8
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, color = 3 
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   'HH + no high order (=ST)', charsize=0.8
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, color = 2
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   'States < 40000 /cm', charsize=0.8			
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.20*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.20*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 0
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.20*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8
   
   
   plot, T, Qmol[i,*], title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'Q' 
   oplot, T, Qmol[i,*], thick = 3
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 5, color = 4
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.05*(!y.crange[1]-!y.crange[0]), $
	   'Irwin (1987)', charsize=0.8
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 2, color = 2
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.10*(!y.crange[1]-!y.crange[0]), $
	   'NextGen', charsize=0.8
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 1, color = 3
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.15*(!y.crange[1]-!y.crange[0]), $
	   'Sauval & Tatum (1984)', charsize=0.8			
   oplot, [!x.crange[0] + 0.60*(!x.crange[1]-!x.crange[0]), $
           !x.crange[0] + 0.68*(!x.crange[1]-!x.crange[0])], $
	  [!y.crange[0] + 0.20*(!y.crange[1]-!y.crange[0]), $
	   !y.crange[0] + 0.20*(!y.crange[1]-!y.crange[0])], $
	   linestyle = 0, thick = 3
   xyouts, !x.crange[0] + 0.70*(!x.crange[1]-!x.crange[0]), $
           !y.crange[0] + 0.20*(!y.crange[1]-!y.crange[0]), $
	   'This work', charsize=0.8	
	      
	   
   ind = where(miNG eq molid[i], nind)
   
   if nind gt 0 then begin
      QNG = 10.d0^poly(alog10(theta), aNG[*,ind])
      oplot, T, QNG, linestyle = 2, color = 2
   endif
   		
   ind = where(miST eq molid[i], nind)
   if nind gt 0 then begin
      QST = 10.d0^poly(alog10(theta), aST[*,ind])
      oplot, T, QST, linestyle = 1, color = 3
      
   endif
   
   ; Irwin 1987
   if molid[i] eq 'H2' then begin
      AI87 = [1.67298118410d4, -1.49945289142d4, 5.74838863349d3, -1.22210505066d3, $
             1.55637569965d2, -1.18744926193d1, 5.02617615447d-1,-9.10563051348d-3]
      QI87 = exp(poly(alog(T), AI87))
      oplot, T, QI87, linestyle = 5, color = 4
   endif
   if molid[i] eq 'CO' then begin
      AI87 = [5.05610415417d4, -5.19600025580d4, 2.33277267148d4, -5.97599449706d3, $
              9.55509531681d2, -9.76517012179d1, 6.22988547018d0, -2.26856284960d-1, $
              3.61025385248d-3]
      QI87 = exp(poly(alog(T), AI87))
      oplot, T, QI87, linestyle = 5, color = 4
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
   
   plot, T, Qmol[i,*], /ylog, /xlog, title = molid[i], xtitle = 'T [K]', $
                        ytitle = 'Q', yr = [ymin, ymax], thick = 3
                        
   oplot, [10, 100], $
	  [10.^(alog10(ymin) + 0.01* alog(ymax/ymin)), 10.^(alog10(ymin) + 0.01* alog(ymax/ymin))], $
	   linestyle = 5, color = 4
   xyouts, 100, $
           10.^(alog10(ymin) + 0.01* alog(ymax/ymin)), $
	   'Irwin (1987)', charsize=0.8
   oplot, [10, 100], $
	  [10.^(alog10(ymin) + 0.04* alog(ymax/ymin)), 10.^(alog10(ymin) + 0.04* alog(ymax/ymin))], $
	   linestyle = 2, color = 2
   xyouts, 100, $
           10.^(alog10(ymin) + 0.04* alog(ymax/ymin)), $
	   'NextGen', charsize=0.8
   oplot, [10, 100], $
	  [10.^(alog10(ymin) + 0.07* alog(ymax/ymin)), 10.^(alog10(ymin) + 0.07* alog(ymax/ymin))], $
	   linestyle = 1, color = 3
   xyouts, 100, $
           10.^(alog10(ymin) + 0.07* alog(ymax/ymin)), $
	   'Sauval & Tatum (1984)', charsize=0.8			
   oplot, [10, 100], $
	  [10.^(alog10(ymin) + 0.10* alog(ymax/ymin)), 10.^(alog10(ymin) + 0.10* alog(ymax/ymin))], $
	   linestyle = 0, thick = 3
   xyouts, 100, $
           10.^(alog10(ymin) + 0.10* alog(ymax/ymin)), $
	   'This work', charsize=0.8
			
			
   ind = where(miNG eq molid[i], nind)
   
   if nind gt 0 then begin
      QNG = 10.d0^poly(alog10(theta), aNG[*,ind])
      oplot, T, QNG, linestyle = 2, color = 2 
   endif	
   
   if molid[i] eq 'H2' then oplot, T, QI87, linestyle = 5, color =4
			
   
   
   if nind2 gt 0 then begin      
      oplot, T, QST, linestyle = 1, color = 3
      
      if total(Qmol[i,*]) gt 0.d0 then $
      plot, T, (Qmol[i,*]-QST)/Qmol[i,*]*100, /xlog, $
         title = molid[i], xtitle = 'T [K]', $
         ytitle = '!4D!3Q [%] (This work - ST)', thick=3, /ys 
   endif else begin
      plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4		 
   endelse
   
   
  
   
   
endfor

device, /close_file
set_plot, 'x'
!p.multi=0

;stop


restore, outdir1 + 'allequil.idl'

set_plot, 'ps'
device, file = 'molecule_equil.ps', ysize=25, yoffset=2
plotthick, 3
!p.charsize=2
!p.font=-1
!p.multi=[0,1,3]

print, ' '
print, 'plotting molecular equilibrium functions in molecule_equil.ps'

for i = 0, nmol-1 do begin 

   ind = where(finite(lgKpmol[i,*]) eq 0, nind)   ; some infinite values in bad cases
   if nind gt 0 then lgKpmol[i,*] = 0.


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
      oplot, T, lgKpNG, linestyle=2, thick=1       
   endif
   
   ind = where(miST eq molid[i], nind)
   if nind gt 0 then begin      
      lgKpST = poly(alog10(theta), bST[*,ind])-theta*DST[ind[0]]
      oplot, T, lgKpST , linestyle=1, thick=7       
   endif
   
   if nind gt 0 then begin
    
   plot, T, (lgKpmol[i,*]- lgKpST), $
            title = molid[i], xtitle = 'T [K]', $
                        ytitle = '!4D!3 log K!dp!n', $
			/xlog, xr=[100, 10000], thick = 5         
   
   plot, T, (lgKpmol[i,*]- lgKpST), $
            title = molid[i], xtitle = 'T [K]', $
                        ytitle = '!4D!3 log K!dp!n', $
			/xlog, xr=[0.1, 10000], thick = 5
   endif else begin
   
      plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4		
   endelse
      		 
   
endfor

device, /close_file
set_plot, 'x'
plotthick, 1
!p.multi=0
 end