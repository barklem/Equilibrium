; 
pro comp_mol_equil, T, molstruct, Qa1, Qa2, Qmol, mu, lgKp, lgKn=lgKn


;compute the equilibrium constant from equation 7 of Tatum 1966
; this means the result is in SI units, N/m^2 = Pa 
; mu should be in amu, molstruct.D in eV, T in K

k = 1.38066e-23 ; SI, J/K

lgKp = alog10(Qa1*Qa2/Qmol) + 2.5d0*alog10(T) + 1.5*alog10(mu) + 3.41405 - 5039.9*molstruct.D/T  ; in partial pressures
Kp = 10.d0^lgKp
lgKn = alog10(Kp/(k*T))   ; in number densities - not tested well and note equation in Tatum seems to be wrong 
; eqn 4 should be I think Kp = kT Kn, which is confirmed by Rossi and Maciel 1983

; compute molar heat capacity
; eqn VIII,15 Herzberg I
;R = 8.31451  ; SI J/mole/K
;Cpo = 2.5*R +  R*deriv(T,T*T*deriv(T,alog(Qtot2)))

end
