pro make_equil

; save, file = outdir1 + 'allpartf.idl', t, Qmol, Qmol4, QmolSt, QmolHH, molid, Dmol

outdir1 = '../molcalc/PartFuncs_Results/'
restore, outdir1 + 'allpartf.idl'

; the atomic data is only valid above 1000 K probably

ind = where(t ge 1000)
t = t(ind)
theta = 5040.d0/T
Qmol = Qmol(*,ind)
QmolST = QmolST(*,ind)


h1 = partitionfunction(t, 1,0, atomfile='atomdata_nist')
c1 = partitionfunction(t, 6,0, atomfile='atomdata_nist')
n1 = partitionfunction(t, 7,0, atomfile='atomdata_nist')
o1 = partitionfunction(t, 8,0, atomfile='atomdata_nist')

h1st = 10.d0^poly(alog10(theta), [0.30103, -.00001])
c1st = 10.d0^poly(alog10(theta), [0.96752, -.09452, 0.08055])
n1st = 10.d0^poly(alog10(theta), [0.60683, -.08674, 0.30565, -0.28114])
o1st = 10.d0^poly(alog10(theta), [0.95033, -.05703])

read_stm, miST, aST, bST, DST

; (CH, NH, OH, C2, CN, CO)

mols = ['H2', 'CH', 'NH', 'OH', 'C2', 'CN', 'CO']
a1 = ['H', 'C', 'N', 'O', 'C', 'C', 'C']
a2 = ['H', 'H', 'H', 'H', 'C', 'N', 'O']
nmols = 7

print, 'T(K)', t, format = '(a30, 50e15.4)'


for i = 0, nmols - 1 do begin
ind = where(molid eq mols[i])
ind2 = where(miST eq mols[i], nind2)
Qmol1 = transpose(Qmol[ind,*])
Qmol1ST = transpose(QmolST[ind,*])
Qmol1STST = 10.d0^poly(alog10(theta), aST[*,ind2])
lgKpST = poly(alog10(theta), bST[*,ind2])-theta*DST[ind2[0]]
Dmol1 = Dmol[ind]
Dmol1 = Dmol1[0]
DST1 = DST(ind2)
DST1 = DST1(0)
if a1[i] eq 'H' then Qa1 = h1
if a1[i] eq 'C' then Qa1 = c1
if a1[i] eq 'N' then Qa1 = n1
if a1[i] eq 'O' then Qa1 = o1
if a2[i] eq 'H' then Qa2 = h1
if a2[i] eq 'C' then Qa2 = c1
if a2[i] eq 'N' then Qa2 = n1
if a2[i] eq 'O' then Qa2 = o1

if a1[i] eq 'H' then Qa1st = h1st
if a1[i] eq 'C' then Qa1st = c1st
if a1[i] eq 'N' then Qa1st = n1st
if a1[i] eq 'O' then Qa1st = o1st
if a2[i] eq 'H' then Qa2st = h1st
if a2[i] eq 'C' then Qa2st = c1st
if a2[i] eq 'N' then Qa2st = n1st
if a2[i] eq 'O' then Qa2st = o1st

if a1[i] eq 'H' then m1 = 1.
if a1[i] eq 'C' then m1 = 12.
if a1[i] eq 'N' then m1 = 14.
if a1[i] eq 'O' then m1 = 16.
if a2[i] eq 'H' then m2 = 1.
if a2[i] eq 'C' then m2 = 12.
if a2[i] eq 'N' then m2 = 14.
if a2[i] eq 'O' then m2 = 16.

mu = (m1 * m2) / (m1 + m2)
comp_mol_equil, T, Dmol1, Qa1, Qa2, Qmol1, mu, lgKp
comp_mol_equil, T, DST1, Qa1st, Qa2st, Qmol1ST, mu, lgKpST1

print, ' '
print, mols(i), 'D', Dmol1, format = '(2a15, 50f15.8)'

print, mols(i), 'D_ST', DST1, format = '(2a15, 50f15.8)'

print, mols(i), 'Q', Qmol1, format = '(2a15, 50e15.4)'
print, mols(i), 'Q_ST', Qmol1ST, format = '(2a15, 50e15.4)'

print, mols(i), 'log10(Kp)', lgKp, format = '(2a15, 50e15.4)'
  
print, mols(i), 'log10(Kp)_ST', lgKpST, format = '(2a15, 50e15.4)'
;print, mols(i), 'log10(Kp)_ST_my', lgKpST1, format = '(2a15, 50e15.4)'

endfor

stop
end