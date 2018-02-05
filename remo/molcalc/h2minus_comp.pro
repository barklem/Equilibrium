; compare h2- partf and equil const

; next gen data
DoNG = 0.7300   
bNG = [19.7162, -5.0018, -2.7680, -1.2845, -0.9859, -0.3380]
aNg = [ 1.0000,  0.0000,  0.0000,  0.0000,  0.0000]

; Sauval & Tatum data
DoST = 0.7300
bST = [11.1759, -0.8735, -0.7470,  0.2748,  0.0000,  0.0000] ;H2 data
aST = [ 1.9508, -1.6265,  0.7472, -0.2751,  0.0000] ; H2- data

T = findgen(100) * 100. + 100
theta = 5040./T
lgKpNG = poly(alog10(theta), bNG)-theta*DoNG
lgKpST = poly(alog10(theta), bST)-theta*DoST

window, 0, retain=2
plot, T, lgKpNG, /xlog, xr=[100, 10000], ytitle = 'log Kp', xtitle = 'T'
oplot, T, lgKpST, psym=4

end
