; 
pro comp_mol_equilkp, T, molstruct, Qa1, Qa2, Qmol, mu, lgKp


;compute the equilibrium constant

lgKp = alog10(Qa1*Qa2/Qmol) + 2.5d0*alog10(T) + 1.5*alog10(mu) + 3.41405 - 5039.9*molstruct.D/T
;Kp = 10.d0^lgKp

; compute molar heat capacity
; eqn VIII,15 Herzberg I
;R = 8.31451  ; SI J/mole/K
;Cpo = 2.5*R +  R*deriv(T,T*T*deriv(T,alog(Qtot2)))

end
