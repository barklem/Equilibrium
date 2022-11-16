; compare h2 partf and equil const

; next gen data
DoNG = 2.6508   
bNG = [15.8052, 33.7578, 34.5956, 27.3455, 16.6214, 9.9717]
aNg = [ 2.5410, -2.4336,  1.4979,  0.0192, -0.7483]

; Sauval & Tatum data
DoST = 2.6508
bST = [ 9.9835, -0.0664, -1.4979, -0.0195,  0.7486] 
aST = [ 2.5410, -2.4336,  1.4979,  0.0192, -0.7483] 

T = findgen(100) * 100. + 100
theta = 5040./T
lgKpNG = poly(alog10(theta), bNG)-theta*DoNG
lgKpST = poly(alog10(theta), bST)-theta*DoST

window, 0, retain=2
plot, T, lgKpNG, /xlog, xr=[100, 10000], yr=[0,80], ytitle = 'log Kp', xtitle = 'T'
oplot, T, lgKpST, psym=4

end
