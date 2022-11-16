pro comp_atom_partf, T, J, E, Qi

; computes partition function for atom
; input J and E (cm-1) for each level
; internal partition func Qi returned for temps T

g = 2*J+1
Ej = E /8065.54 * 1.602d-19  ; cm-1 > J
k = 1.38066d-23

Qi = 0.d0 * T

for i = 0, n_elements(J) -1 do begin
  Qi_c = g[i]*exp(-Ej[i]/k/T)
  Qi = Qi + Qi_c
  print, i, J[i], g[i], E[i], Ej[i], -Ej[i]/k/T, Qi_c, Qi_c/Qi, Qi
endfor
print, 'Qi = ', Qi
stop
return
end
