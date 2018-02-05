function WCALC, XNH, XNE, XNHE, NS, TEMP
;  	
;  Computes the occupation probability for a H atom in 
;  state with effective principle quantum number NS in a plasma
;  enviroment with NH, NE, NHE number densities of H, ions,
;  and He atoms respectively.  This code assumes the perturbing
;  neutral H and He are in the ground state, (noting the hard
;  sphere approximation used is quite crude anyway) and that ions 
;  are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
;
;  See eqn 4.71 Hummer & Milhalas (1988) 
;  Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
;  which is a fit to detailed calculation including correlations,
;  thus the temperature dependence
;
;  Sizes of neutral atoms adopted are sqrt(<r^2>)
;  
;  Coded by Paul Barklem and Kjell Eriksson, Aug 2003
; 	
     
      IONH=2.17991E-11
      A0=5.29177E-9
      E=4.803207E-10
;
      NS2 = NS*NS
      NS4 = NS2*NS2
;  
;  Neutral perturbers
;
      CHI=IONH/(NS2)
      RIH=SQRT(2.5*NS4 + 0.5*NS2)*A0
      X1=RIH + 1.73*A0
      X2=RIH + 1.02*A0
      NEUTR=XNH*X1*X1*X1 + XNHE*X2*X2*X2
      WNEUTR = EXP(-4.18879 * NEUTR)  ; 4.18879 is 4*!pi/3
;
;  Charged perturbers
;      
      K=1.
      IF (NS gt 3.) THEN BEGIN
;       K=5.33333333*(NS/(NS+1.))^2 * (NS + 7./6.)/(NS2+NS+0.5)
        K=5.33333333*NS/(NS+1.)/(NS+1.)
      ENDIF
;      ION  = NE*16.*(E*E/CHI/DSQRT(K)^3
;      WION = DEXP(-4.*!PI/3. * (ION))
      IF ((XNE gt 10.) AND (TEMP gt 10.)) THEN BEGIN  ; just in case!!! 
        A = 0.09 * EXP(0.16667d0*ALOG(XNE)) / SQRT(TEMP)
        X = EXP(3.15*ALOG(1.+A))
        BETAC = 8.3E14 * EXP(-0.66667*ALOG(XNE)) * K / NS4
        F = 0.1402d0*X*BETAC*BETAC*BETAC /(1.+0.1285*X*BETAC*SQRT(BETAC))
        WION = F/(1.d0+F)
      ENDIF ELSE BEGIN 
        WION = 1.0d0     
      ENDELSE
;
      WCALC = WION * WNEUTR 
      RETURN, WCALC
      END
