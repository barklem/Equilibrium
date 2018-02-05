pro plotcut_exp



xns = dindgen(299)/10. + .1
n = n_elements(xns)
w1 = xns*0.
w1b = xns*0.
w2 = xns*0.
w3 = xns*0.

for i = 0, n-1 do begin
   w1[i] = wcalc(1e14,1e10,1e13,xns[i],4000.)
   w1b[i] = wcalc(4e16,4e12,4e15,xns[i],5000.)
   w2[i] = wcalc(1e17,3e13,1e16,xns[i],6000.)
   w3[i] = wcalc(1e17,7e14,1e16,xns[i],8000.)
endfor


; plot parameters ps files
encapsulated  = 1
landscape     = 0
font          = -1
default_thick = 4
charsize      = 1.2
charthick     = default_thick


color_1    = 'blu4'
color_2   = 'tomato'
color_3    = 'grn4'
color_4   = 'plum'

cgps_close

cgPS_open, 'wplot_exp.eps', encapsulated=encapsulated, landscape=landscape, $
                 font=font, charsize=charsize, default_thick=default_thick, $
                 /nomatch
cgplot, xns, w2, col=color_1, thick=7, ytitle = 'Occupation probability', xtitle = 'effective principal quantum number n!u*!n', xr=[0,8], yr=[0.98,1.002], /ys, /xs
cgplot, /ov, xns, w1, linest = 5, col=color_2, thick=7
cgplot, /ov, xns, w1b, linest = 3, col=color_4, thick=7
cgplot, /ov, xns, w3, linest = 2, col=color_3, thick=7


xloc   = !x.crange[0]+(!x.crange[1]-!x.crange[0])*0.1
yloc   = !y.crange[0]+(!y.crange[1]-!y.crange[0])*0.3
cgLegend, Title=['4000 K', '5000 K', '6000 K', '8000 K'], $
                Lines=[5,3,0,2], thick=7, Color=[color_2,color_4,color_1,color_3], $
                charsize=charsize, charthick=charthick, vspace=1.8, length=0.125, $
                Location=[xloc,yloc], /Data, /background

cgps_close


; convert EPS files to PDF
spawn,'for i in *'+'.eps'+'; do epstopdf $i ; done'

end
