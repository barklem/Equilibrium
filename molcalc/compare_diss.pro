read_diss, name, adopted, hh=hh, luo=luo, g2=g2, other=other

;n=n_elements(adopted)

ps = 1
if ps eq 0 then begin
   set_plot, 'x'
   window, 0, ysize = 1000, xsize = 1500
   !p.charsize = 3
endif else begin
   set_plot, 'ps'
   device, file = 'Dissociation/diss.eps', /landscape
   !p.charsize = 1.5
endelse

!p.multi = [0,3,2]

psymcircle
ind = where(luo(*,0) gt 0.01 and hh(*,0) gt 0.01, nind)
plot, luo(ind,0), hh(ind,0), psym = 8, ytitle = 'Huber & Herzberg (1979) - compilation', xtitle = 'Luo (2007) - Expt. compilation'
oplot, [0,20], [0,20], linestyle = 1
for i = 0, nind-1 do begin
   if luo(ind(i),1) ne 0.d0 then begin
       oplot, [luo(ind(i),0)-luo(ind(i),1), luo(ind(i),0)+luo(ind(i),1)], [hh(ind(i),0), hh(ind(i),0)]
   endif
endfor
for i = 0, nind-1 do begin
   if hh(ind(i),1) ne 0.d0 then begin
       oplot, [luo(ind(i),0), luo(ind(i),0)], [hh(ind(i),0)-hh(ind(i),1), hh(ind(i),0)+hh(ind(i),1)]
   endif
endfor

ind = where(luo(*,0) gt 0.01 and g2(*,0) gt 0.01, nind)
plot, luo(ind,0), g2(ind,0), psym = 8, ytitle = 'Curtiss et al (1991) - Gaussian-2 theory', xtitle = 'Luo (2007) - Expt. compilation' 
oplot, [0,20], [0,20], linestyle = 1
for i = 0, nind-1 do begin
   if luo(ind(i),1) ne 0.d0 then begin
       oplot, [luo(ind(i),0)-luo(ind(i),1), luo(ind(i),0)+luo(ind(i),1)], [g2(ind(i),0), g2(ind(i),0)]
   endif
endfor

ind = where(g2(*,0) gt 0.01 and hh(*,0) gt 0.01, nind)
plot, g2(ind,0), hh(ind,0), psym = 8, xtitle = 'Curtiss et al (1991) - Gaussian-2 theory', ytitle = 'Huber & Herzberg (1979) - compilation'
oplot, [0,20], [0,20], linestyle = 1
for i = 0, nind-1 do begin
   if hh(ind(i),1) ne 0.d0 then begin
       oplot, [g2(ind(i),0), g2(ind(i),0)], [hh(ind(i),0)-hh(ind(i),1), hh(ind(i),0)+hh(ind(i),1)]
   endif
endfor

print, 'Luo vs HH'
print, '        Luo    HH     diff'
diff = hh(*,0) - luo(*,0) 
diff1 = abs(diff)
ind2 = where(hh(*,0) ne 0.d0 and luo(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
print, name(ind2(ind(i))), luo(ind2(ind(i)),0), hh(ind2(ind(i)),0), diff(ind2(ind(i))), diff(ind2(ind(i)))/luo(ind2(ind(i)),0), format = '(a5, 4f7.2)'
endfor

plot, luo(ind2,0), diff(ind2)/luo(ind2,0), psym = 8, /ys, yr=[-0.5, 0.5], xtitle = 'Luo (2007)', ytitle = 'relative difference'

print, 'Luo vs G2'
print, '        Luo    G2     diff'
diff =  g2(*,0) - luo(*,0) 
diff1 = abs(diff)
ind2 = where(g2(*,0) ne 0.d0 and luo(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
print, name(ind2(ind(i))), luo(ind2(ind(i)),0), g2(ind2(ind(i)),0), diff(ind2(ind(i))),  diff(ind2(ind(i)))/luo(ind2(ind(i)),0), format = '(a5, 4f7.2)'
endfor


plot, luo(ind2,0), diff(ind2)/luo(ind2,0), psym = 8, /ys, yr=[-0.5, 0.5], xtitle = 'Luo (2007)', ytitle = 'relative difference'

print, 'HH vs G2'
print, '        HH     G2     diff'
diff = hh(*,0) - g2(*,0)
diff1 = abs(diff)
ind2 = where(g2(*,0) ne 0.d0 and hh(*,0) ne 0.d0, nind2)
ind = reverse(sort(diff1(ind2))) 
for i = 0, nind2-1 do begin
print, name(ind2(ind(i))), hh(ind2(ind(i)),0), g2(ind2(ind(i)),0), diff(ind2(ind(i))),  diff(ind2(ind(i)))/luo(ind2(ind(i)),0), format = '(a5, 4f7.2)'
endfor


plot, g2(ind2,0), diff(ind2)/g2(ind2,0), psym = 8, /ys, yr=[-0.5, 0.5], xtitle = 'Curtiss et al (1994)', ytitle = 'relative difference'

if ps ne 0 then begin
   device, /close_file
   set_plot, 'x'
endif


end