tio_np = double([6.8700, 11.9229, -1.4044, 0.7899, -0.7317, -0.0193, -0.4994])
tio_st = double([6.8700, 11.4047, -1.1484, 0.6478, -0.6737])

temp = dindgen(1001)/1000*5900+100.
th = 5040.d0/temp
tio_kp_np = poly(alog10(th), tio_np[1:6])-th*tio_np[0]
tio_kp_st = poly(alog10(th), tio_st[1:4])-th*tio_st[0]

; NIST
t = [100., 200., 300., 500., 600., 1000., 2000., 3000., 4000., 5000., 6000.]
kf = 10.^[-23.423, -9.008, -4.240, -0.497, 0.415, 2.173, 3.219, 3.255, 2.613, 1.439, 0.630]

; tatum66
t2 = [1000., 2000., 3000., 4000., 5000., 6000.]
part2 = [1.027e+4, 3.081e+4, 6.308e+4, 1.087e+5, 1.707e+5, 2.531e+5]
; atomic partition functions from Sauval & Tatum
; probably bad at low T
lgth = alog10(5040.d0/t2)
Qoxy = 10.^(0.95033 - 0.05703*lgth)
Qtit = 10.^(1.47343 - 0.97220*lgth + 1.47986*lgth*lgth - 0.93275*lgth*lgth*lgth)
;compute the equilibrium constant
D = 6.8
moxy = 15.9994  ; oxygen mass in amu
mtit = 47.90    ; titanium mass in amu
mu = (moxy*mtit)/(moxy+mtit)    ; reduced mass of molecule in amu
lgKp = alog10(Qoxy*Qtit/Qtot) + 2.5*alog10(t2) + 1.5*alog10(mu) + 3.41405 - 5039.9*D/t2
kp2 = lgKp

k = 1.38066d-23
R = 8.314472d-3  ; kJ/mole K

window, 0, retain=2
plot, temp, tio_kp_np, xr=[10, 6100], /xs
oplot, temp, tio_kp_st, linestyle = 2
;oplot, t, kp, psym=5
;oplot, t, alog10(kf*exp(-6.87*1.602e-19/k/t)*k*t), psym=4
oplot, t, alog10((1./kf)*exp(-6.87*1.602e-19/k/t)/k/t/1e5), psym=2, thick=5
;oplot, t, -alog10(kf)-(6.87*1.602e-19/k/t)*alog10(2.7)-alog10(k*t)-5., psym=2, thick=5
;oplot, t, -alog10(kf)-5034.*6.87/t-alog10(t)+17.86, psym=2, thick=5

oplot, t2, kp2, psym=2
stop
end
