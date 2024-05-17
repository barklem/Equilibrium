; 
pro comp_mol_equil3, T, molstruct, Qa, Qb, Qmol, mu, Ib, lgKp, lgKn=lgKn

; compute the 3-body equilibrium constant for a negative molecular ion
; i.e. A + B + e <-> AB-
; assuming the dissocation energy molstruct.D is for
; AB- -> A + B-, then the ionisation energy Ib should be for B-
;  i.e. the electron affinity of B

; the result is in SI units, N/m^2 = Pa 
; mu should be in amu, molstruct.D and Ib in eV, T in K

k = 1.38066e-23 ; SI, J/K

lgKp = alog10(2.*Qa*Qb/Qmol) + 5d0*alog10(T) + 1.5*alog10(mu) + 1.9371 - 5039.9*(molstruct.D+Ib)/T  ; in partial pressures
Kp = 10.d0^lgKp
lgKn = alog10(Kp/((k*T)^2.d0))   ; in number densities 

end
