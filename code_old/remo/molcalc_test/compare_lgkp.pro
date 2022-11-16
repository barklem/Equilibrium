pro compare_lgkp, label,co=co, debug=debug


me	= 9.109382e-28
kk	= 1.3806504e-16
hh	= 6.62606896e-27
evtoerg	= 1.602176487e-12
amu	= 1.660538782e-24                       


restore,'~/People/barklem/molcalc_test/alldata.sav'
ww=where(t ge 1000. and t le 16000.)
adiv=20.0 & xt=alog10(t)-4.1 & div=(sqrt(1.+(adiv*xT/2.)^2)+adiv/2.*xT)/adiv*5.0
th=5039.7475/t

molid_m=['H2', 'H2+', 'H2O', 'OH', 'CH',  'CO', 'CN', 'C2', 'N2', 'O2', 'NO', 'NH']

aka=[12.5277d0,10.3642d0,25.4200d0,12.2799d0,11.6371d0,13.9708d0,12.8372d0,12.7523d0,14.0287d0,13.0168d0,13.0100d0,11.8225d0]
ak1=[4.72813d0,0.980461d0,10.5220d0,5.02288d0,3.35791d0,12.3989d0,8.45179d0,6.73465d0,12.1756d0,5.13878d0,7.96834d0,3.81278d0]
ak2=[-1.92976d-1,-1.53083d0,1.69390d-1,1.70121d-1,-3.37335d-1,8.01480d-1,5.00672d-1,1.91824d-1,1.88811d0,-2.21281d-1, 9.90928d-1,1.79164d-2]
ak3=[-1.09505d-1,-5.56921d-1,1.8368d-2,3.40739d-2,-1.40312d-1,2.77516d-1,2.09641d-1,5.58272d-2,7.23205d-1,-9.79568d-2,3.54155d-1,-1.71642d-2]
ak4=[-1.74118d-2,-7.51944d-2,8.1730d-4,3.06532d-3,-2.01630d-2,-1.74118d-2,-7.51944d-2,8.1730d-4,3.06532d-3,-2.01630d-2,4.83538d-2,-3.69515d-3]




read_stm, molid_st, aks, bks, dks



if (label eq 'H2O') then begin
 a_oi=[-6.26219795D+00,  5.12651680D+00, -1.44691306D+00,  2.23196759D-01, -1.78196943D-02,  5.71458469D-04]                   
 a_hi=[-2.61655891D+02,  1.63428326D+02, -4.06133526D+01,  5.03282928D+00, -3.10998364D-01,  7.66654594D-03]                   
 a_h2o=[1.44518994D+03, -1.18568100D+03,  4.02064308D+02, -7.20425879D+01,  7.19907707D+00, -3.80370337D-01, 8.31134943D-03]   
 lnq_oi=poly(alog(t),a_oi)
 lnq_hi=poly(alog(t),a_hi)
 lnq_h2o=poly(alog(t),a_h2o)
 d_h2o=9.512       


 const=(2.*!DPI*amu/hh/hh)                                                                          

 lnkp_h2o=2.*lnq_hi+lnq_oi-lnq_h2o+ 3.*alog(const)+5.*alog(kk*t)-(d_h2o*evtoerg/kk/t)

 m   = where(molid_m  eq label)
 akd= (aka(m[0])-(ak1(m[0])-(ak2(m[0])-(ak3(m[0])-ak4(m[0])*th)*th)*th)*th)
 
 plot,t[ww],akd[ww],th=2,xtit='Temp',ytit='Log10(Kp)',tit=label
 oplot,t[ww],lnkp_h2o[ww]/alog(10.),th=2,co=2000

 co=poly_fit(th[ww],lnkp_h2o[ww]/alog(10.),4)

endif else begin

;
; compute and plot dissoc constants
;
 wid = where(molid    eq label)
 m   = where(molid_m  eq label)
 i   = where(molid_st eq label)

 akd= (aka(m[0])-(ak1(m[0])-(ak2(m[0])-(ak3(m[0])-ak4(m[0])*th)*th)*th)*th)
 akds= poly(alog10(th),bks[*,i[0]])-dks[i[0]]*th

 plot,t[ww],akd[ww],th=2,xtit='Temp',ytit='Log10(Kp)',tit=label
 oplot,t[ww],lgkpmol[wid[0],ww],co=2000,th=2
 oplot,t[ww],akds[ww],ps=8,co=7000

 co=poly_fit(th[ww],lgkpmol[wid[0],ww],4)

 co[0]=co[0]+1.					; SI --> cgs
 
endelse

print, 'POLY_FIT LOG10(KP(TH)): ',co

if keyword_set(debug) then stop
end
