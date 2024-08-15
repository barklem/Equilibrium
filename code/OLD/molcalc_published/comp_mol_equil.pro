; 
pro comp_mol_equil, T, molstruct, Qa1, Qa2, Qmol, mu, lgKp, lgKn=lgKn


;compute the equilibrium constant

k = 1.38066e-16 ; cgs

lgKp = alog10(Qa1*Qa2/Qmol) + 2.5d0*alog10(T) + 1.5*alog10(mu) + 3.41405 - 5039.9*molstruct.D/T  ; in partial pressures
Kp = 10.d0^lgKp
lgKn = alog10(Kp*k*T)   ; in number densities

; compute molar heat capacity
; eqn VIII,15 Herzberg I
;R = 8.31451  ; SI J/mole/K
;Cpo = 2.5*R +  R*deriv(T,T*T*deriv(T,alog(Qtot2)))

end
