; calculate and compare atomic partition functions

filenm = 'ti_nist.dat'
;filenm = 'o_nist.dat'
filenm = 'c_nist.dat'
;filenm = 'n_nist.dat'

T = [1., 2., 5., 10., 50., 100., 200., 500., 1000., 2000., 5000., 8000.]

comp_atompart, filenm, T, Qi 

print, T, format='(15i10)'
print, Qi, format='(15e10.3)'

tt = findgen(1000)*10.+10.
lgth = alog10(5040.d0/tt)
Qoxy = 10.^(0.95033 - 0.05703*lgth)
Qtit = 10.^(1.47343 - 0.97220*lgth + 1.47986*lgth*lgth - 0.93275*lgth*lgth*lgth)
Qcar = 10.^(0.96752 - 0.09452*lgth + 0.08055*lgth*lgth)
Qnit = 10.^(0.60683 - 0.08674*lgth + 0.30565*lgth*lgth - 0.28114*lgth*lgth*lgth)

psymcircle
plot, T, Qi, psym=8, /ylog, /xlog, xr=[0.5, 10000]
oplot, tt, Qcar

stop

end


 
