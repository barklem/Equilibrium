read_mol, 'H2_nist_huber79.mol', molstruct
T = 10.^(findgen(501)/100-1.)
comp_mol_partf, T, molstruct, Qi
comp_mol_partf_homonuc, T, molstruct, 0.5, Qi2

Bv = molstruct.Be[0] - molstruct.alfe[0]*0.5 + molstruct.game[0]*0.25  
splt = Bv*2
h = 6.626076d-27   ; cgs
c = 2.997924d10
k = 1.38066d-16
hckT = h*c/k/T
splt_e = exp(-hckT*splt)
splt_e3 = exp(-hckT*splt*2)
splt_e4 = exp(-hckT*splt*3)
splt_e5 = exp(-hckT*splt*4)

pop1 = 0.25 / Qi2
pop2 = 0.75 * 3. * splt_e / Qi2
pop3 = 0.25 * 5. * splt_e3 / Qi2
pop4 = 0.75 * 7. * splt_e4 / Qi2
pop5 = 0.25 * 9. * splt_e5 / Qi2
popt = pop1 + pop2 + pop3 + pop4 + pop5

pop1d = 1. / (Qi*2)
pop2d = 3. * splt_e / (Qi*2)
pop3d = 5. * splt_e3 / (Qi*2)
pop4d = 7. * splt_e4 / (Qi*2)
poptd = pop1d + pop2d + pop3d + pop4d

set_plot, 'ps'
device, file = 'h2_tests.ps'
plot, T, Qi, /ylog, /xlog
oplot, T, Qi2, linestyle=1, thick =3
oplot, T, Qi2*4, linestyle=1, thick =1
plot, T, pop1, /xlog, yr = [1e-10, 2], /ys
oplot, T, pop2
oplot, T, pop3
oplot, T, pop4
oplot, T, popt, thick = 5
oplot, T, pop1d, linestyle = 1, thick = 3
oplot, T, pop2d, linestyle = 1, thick = 3
oplot, T, pop3d, linestyle = 1, thick = 3
oplot, T, pop4d, linestyle = 1, thick = 3
oplot, T, poptd, linestyle = 1, thick = 9
plot, T, poptd / popt, /xlog
device, /close_file


set_plot, 'x'