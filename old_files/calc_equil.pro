; read data

;filenm = 'tio_tatum66.dat'
;filenm = 'tio_nist.dat'
;filenm = 'tio_nist_huber79.dat'
;filenm1 = 'o_nist.dat'
;filenm2 = 'ti_nist.dat'
;m1 = 15.9994  ; oxygen mass in amu
;m2 = 47.90    ; titanium mass in amu
;st = [5.3051, -2.3739, 0.8940, -0.3641]  ; tio partf from sauval and tatum
;stk = [6.8700, 11.4047, -1.1484, 0.6478, -0.6737] ; Kp from sauval and tatum
;outfile = 'TiO.sav'


filenm = 'h2+_nist_huber79.dat'
filenm1 = ' '
filenm2 = ' '
m1 = 1.008  ; H mass in amu
m2 = 1.008  ; H mass in amu
st = [2.5410, -2.4336, 1.4979, 0.0192, -0.7483]  ; h2+ partf from sauval and tatum
stk = [2.6508, 9.9835, -0.0664, -1.4979, -0.0195, 0.7486] ; Kp from sauval and tatum
outfile = 'H2+.sav'
plotfile = 'H2+.eps'

;filenm = 'co_nist_huber79.dat'
;filenm1 = 'o_nist.dat'
;filenm2 = 'c_nist.dat'
;m1 = 15.9994  ; oxygen mass in amu
;m2 = 12.011   ; carbon mass in amu
;st = [3.6076, -1.7608, 0.4172]  ; co partf from sauval and tatum
;stk = [11.092, 12.2263, -0.8829, -0.1230, -0.3226] ; Kp from sauval and tatum
;outfile = 'CO.sav'


;filenm = 'cn_nist_huber79.dat'
;filenm1 = 'c_nist.dat'
;filenm2 = 'n_nist.dat'
;m1 = 12.011   ; carbon mass in amu
;m2 = 14.007   ; nitrogen mass in amu
;st = [4.0078, -2.1514, 0.9226, -0.1671]  ; cn partf from sauval and tatum
;stk = [7.760, 11.4479, -0.4840, -0.4160, -0.9435, 0.8380] ; Kp from sauval and tatum
;outfile = 'CN.sav'

print, filenm
print, filenm1
print, filenm2


openr, 1, filenm
D = 0.
sigma = 0
ns = 0
s = ' '
readf, 1, s ; comment line
readf, 1, D, sigma, ns
readf, 1, s ; comment line

;ns = 1

state = strarr(ns)
Te = fltarr(ns)
gel = intarr(ns)
Be = fltarr(ns)
alfe = fltarr(ns)
we = fltarr(ns)
wxe = fltarr(ns)

for i = 0, ns - 1 do begin
  s = ' '
  readf, 1, s, format='(a100)'
  if strmid(s,0,1) eq '#' then begin
    i = i - 1
    goto, comment
  endif
  s = strsplit(s,' ', /extract)
  state[i] = s[0]
  Te[i] = float(s[1])
  gel[i] = fix(s[2])
  Be[i] = float(s[3])
  alfe[i] = float(s[4])
  we[i] = float(s[5])
  wxe[i] = float(s[6])
;  print, Te[i]
  comment:
endfor

close, 1

; define constants

h = 6.626076d-27   ; cgs
c = 2.997924d10
k = 1.38066d-16

T = 5000.
T = [0.01, 0.1, 1., 3., 5., 10., 30., 50., 100., 300., 500., 1000., 3000., 5000., 10000.]
T = 10.^(findgen(101) * 0.06 - 2.)

hckT = h*c/k/T

Gzero = (0.5*we[0] - 0.25*wxe[0])   ;vibrational energy of lowest state in wavenumbers
fac1 = exp(hckT * Gzero)

Qtot = 0.d0
Qtot2 = 0.d0
Qe = 0.d0

for i = 0, ns-1 do begin
   
   vmax = fix(we[i]/wxe[i]/2. - 0.5)
   
   for v = 0, vmax do begin
   
      ; precompute some oft used terms 
      vph = v + 0.5
      Bv = (Be[i] - alfe[i]*vph)
      
      ; the electronic and vibrational parts in the summation
      ; putting fac1 inside here helps avoid machine limits
      elvib_part = gel[i] * exp( -hckT * (we[i]*vph - wxe[i]*vph*vph + Te[i]) +hckT * Gzero)
      
      ; calculate rotational part, in high T approximation
      Qrot = 1.d0 / sigma / hckT / Bv
      
      ; calculate rotational part, explicitly
      Qrot2 = dblarr(n_elements(T))*0.
      for k = 0, n_elements(T)-1 do begin
        Jmax = sqrt(1./(2*Bv*hckT[k]))    ; J where function is maximum
	Jmax = max([Jmax, 10])              ; very low T, Jmax -> 0
        J = findgen(fix(Jmax*10))         ; 10 times this maximum seems to suffice
        Qrot2[k] = total((2*J+1)*exp(-hckT[k]*Bv*J*(J+1))) / sigma
      endfor
      
      ; here's an approximate vibrational partition function
;      Qvib = 1./(1.-
          
      Qtot = Qtot + Qrot * elvib_part
      Qtot2 = Qtot2 + Qrot2 * elvib_part
      
;      print, T, Qrot2/Qrot
;       print, v, -hckT * (we[i]*vph - wxe[i]*vph*vph + Te[i]),  -hckT 
      
;      print, v, Qtot, Qrot, elvib_part

   endfor
   
;   print, '**', i, Qtot
endfor

;Qtot = fac1 * Qtot
;Qtot2 = fac1 * Qtot2

; atomic partition functions 

if filenm1 ne ' ' then begin
comp_atompart, filenm1, T, Q1 
endif else begin
Q1 = 10.^(0.30103 - 0.00001*alog10(5040.d0/T))
endelse

if filenm2 ne ' ' then begin
comp_atompart, filenm2, T, Q2 
endif else begin
Q2 = 1.
endelse

;compute the equilibrium constant

mu = (m1*m2)/(m2+m1)    ; reduced mass of molecule in amu
lgKp = alog10(Q1*Q2/Qtot) + 2.5*alog10(T) + 1.5*alog10(mu) + 3.41405 - 5039.9*D/T
lgKp2 = alog10(Q1*Q2/Qtot2) + 2.5*alog10(T) + 1.5*alog10(mu) + 3.41405 - 5039.9*D/T
Kp = 10.^lgKp
Kp2 = 10.^lgKp2

; compute molar heat capacity
; eqn VIII,15 Herzberg I
R = 8.31451  ; SI J/mole/K
Cpo = 2.5*R +  R*deriv(T,T*T*deriv(T,alog(Qtot2)))


th = 5040.d0/T
Qst_a = 10.^poly(alog10(th), st)
nn = n_elements(stk) - 1
lgKpst_a = poly(alog10(th), stk[1:nn])-th*stk[0]



;print, 'T[K]', T, format = '(a10,15f13.2)'
;print, ' '
;;print, 'Qtot', Qtot, format = '(a10,15e13.5)'
;print, 'Qint', Qtot2, format = '(a10,15e13.5)'
;print, 'Qint(ST)', Qst_a, format = '(a10,15e13.5)'
;print, 'rel err%', (-Qtot2 + Qst_a)*100./Qtot2, format = '(a10,15f13.2)'
;print, ' '
;;print, 'Kp', Kp, format = '(a10,15e10.2)'
;;print, 'Kp2', Kp2, format = '(a10,15e10.2)'
;;print, 'lgKp', lgKp, format = '(a10,15e10.2)'
;print, 'logKp', lgKp2, format = '(a10,15e13.5)'
;print, 'logKp(ST)', lgKpst_a, format = '(a10,15e13.5)'
;print, 'diff', (lgKp2 - lgKpst_a), format = '(a10,15e13.5)'
;;print, 'Cpo2', Cpo, format = '(a10,15e10.2)'
;;print, 'Cpo(int)', Cpo-2.5*R, format = '(a10,15e10.2)'

;print, 'T[K]', 




Tst = findgen(10000)+1.
th = 5040.d0/Tst
Qst = 10.^poly(alog10(th), st)
nn = n_elements(stk) - 1
lgKpst = poly(alog10(th), stk[1:nn])-th*stk[0]

set_plot, 'ps'
device, file = plotfile
;!p.multi=[0,1,2]

psymcircle
plot, T, Qtot2, psym=8, /ylog, /xlog, xr=[0.005,12000.], /xs, $
      ytitle = 'Q!dint!n', xtitle = 'T[K]'
;oplot, T, Qtot, psym=6
oplot, T, Qst_a

plot, T, lgKp2, psym=8, /xlog, xr=[100.,10000.], /xs, $
      ytitle = 'log K!dp!n', xtitle = 'T[K]'
oplot, T, lgKpst_a

device, /close_file
set_plot, 'x'


Q = Qtot2
Q_st = Qst_a
lgKp = lgKp2
lgKp_st = lgKpst_a

save, filename=outfile, T, Q, lgKp, Q_st, lgKp_st

stop
end
