; try and compute equilibrium constant for H as a check

t = [100., 1000., 5000.]
I = 13.6
u1 = 1.
u0 = 2.
logK = -5040./t*i + 2.5*alog10(t) + alog10(u1/u0) - 0.1762
K2 = (2.d0*!pi*9.10953e-31)^1.5*(1.38e-23*t)^2.5 / (6.626e-34^3) * 2*u1/u0*exp(-13.6*1.6e-19/1.38e-23/t)
print, t
print, logK
print, alog10(K2)
end
