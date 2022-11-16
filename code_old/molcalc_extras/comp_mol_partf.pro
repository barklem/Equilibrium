;
pro comp_mol_partf, T, molstruct, Qi, Qi4=Qi4, QiST=QiST, errorfile=errorfile, showinfo=showinfo, Ezero=Ezero

; Qi4 is the partition function only including states below 40000 cm^-1
; 
; QiST retains only terms considered by Sauval and Tatum including only states below 40000 cm^-1
; i.e. wye = 0, game = 0, De = 0, bete = 0
;
;
;  Ezero returns the energy of the lowest state wrt the bottom of potential well
;

ms = molstruct ; for brevity, but retaining clarity in input line
               ; also allows us to temporarily change values

h = 6.626076d-27   ; cgs
c = 2.997924d10
k = 1.38066d-16

hckT = h*c/k/T

Qi = 0.d0
Qi4 = 0.d0
QiST = 0.d0

Qtot = 0.d0
Qe = 0.d0
Qtot4 = 0.d0
QtotST = 0.d0

grnd = 0

if ms.ns le 0 then goto, skipall

; fill some data
;ind = where (ms.we eq 0., nind, complement=ind2, ncomplement=nind2)
;if nind gt 0 and nind2 gt 0 then ms.we[ind] = mean(ms.we[ind2]) 
;ind = where (ms.Be eq 0., nind, complement=ind2, ncomplement=nind2)
;if nind gt 0 and nind2 gt 0 then ms.Be[ind] = mean(ms.Be[ind2]) 


ilow = where(ms.Te eq 0., nlow)
if nlow gt 0 then begin  ; ilow must be a scalar, not an array to work properly below, otherwise Gzero is an array
    Jzero = ms.lambda[ilow]  ; find lowest state
    Bv_zero = ms.Be[ilow] - ms.alfe[ilow] * 0.5 + ms.game[ilow] * 0.25
    Dv_zero = ms.De[ilow] - ms.bete[ilow] * 0.5
    Fzero = Bv_zero * Jzero*(Jzero + 1) - Dv_zero* (Jzero*(Jzero + 1)) * (Jzero*(Jzero + 1))
    qqq = min(Fzero, ilow2)
    ilow = ilow[ilow2]
    if ilow ne ms.ns-1 and keyword_set(errorfile) then begin
           openw, lun, errorfile, /get_lun, /append
           printf, lun, ms.molid + ':   ilow is not last state - ilow = '+string(ilow, '(i2)') +' ns-1 = ' +string(ms.ns-1)
           close, lun
           free_lun, lun
    endif 
endif else begin
  ilow = ms.ns-1   ; this allows code to run anyway, but will flag an error to logfile below if requested
endelse
;print, ilow, ms.we[ilow], ms.wxe[ilow], ms.wye[ilow]
Gzero = (0.5*ms.we[ilow] - 0.25*ms.wxe[ilow] + 0.125*ms.wye[ilow])   ;vibrational energy of lowest state in wavenumbers
GzeroST = (0.5*ms.we[ilow] - 0.25*ms.wxe[ilow])   

; calculate lowest rotational energy. 
Fzero = 0.d0
FzeroST = 0.d0
Fzero2 = 0.d0
Fzero2ST = 0.d0
Jzero = ms.lambda[ilow]
Jzero2 = ms.lambda[ilow] + 1
if Jzero gt 0 then begin
   Bv_zero = ms.Be[ilow] - ms.alfe[ilow] * 0.5 + ms.game[ilow] * 0.25
   Bv_zeroST = ms.Be[ilow] - ms.alfe[ilow] * 0.5 
   Dv_zero = ms.De[ilow] - ms.bete[ilow] * 0.5
   Fzero = Bv_zero * Jzero*(Jzero + 1) - Dv_zero* (Jzero*(Jzero + 1)) * (Jzero*(Jzero + 1))  ; rotational energy of lowest state
   FzeroST = Bv_zeroST* Jzero*(Jzero + 1) 
endif

 

Ezero = Gzero + Fzero 
EzeroST = GzeroST + FzeroST

; we have to adjust if homonuclear, I=0 and lowest rot state doesn't exist
if ms.lambda[ilow] eq 0 and ms.homonuc eq 1 and ms.I1 eq 0 then begin
    
    if ((ms.plusmin[ilow] eq '-' and ms.oddeven[ilow] eq 'g') or (ms.plusmin[ilow] eq '+' and ms.oddeven[ilow] eq 'u')) then begin

         Bv_zero = ms.Be[ilow] - ms.alfe[ilow] * 0.5 + ms.game[ilow] * 0.25
         Bv_zeroST = ms.Be[ilow] - ms.alfe[ilow] * 0.5 
         Dv_zero = ms.De[ilow] - ms.bete[ilow] * 0.5
         Fzero2 = Bv_zero * Jzero2*(Jzero2 + 1) - Dv_zero* (Jzero2*(Jzero2 + 1))  * (Jzero2*(Jzero2 + 1))  ; rotational energy of second lowest state
         Fzero2ST = Bv_zeroST* Jzero2*(Jzero2 + 1)
         Ezero  = Gzero + Fzero2
         EzeroST = GzeroST + Fzero2ST
endif
endif

;if ms.molid eq 'AlH+' then stop

;print, nlow, ilow, ms.we[ilow], ms.wxe[ilow], ms.wye[ilow], Gzero

for i = 0, ms.ns-1 do begin

    if ms.smult[i] lt 0 then ms.smult[i] = 1   ; we assume 1Sigma state if no other info
    if ms.lambda[i] lt 0 then ms.lambda[i] = 0
    
    if ms.Te[i] lt 0. or ms.smult[i] lt 0 or ms.lambda[i] lt 0 then goto, skip  ; no energy level or symmetry
    if ms.we[i] eq 0. then goto, skip  ; no info on vibrational states - note if no mean estimate, then all zero and should lead to skipping molecule
    if ms.Te[i] eq 0. and ms.homonuc ne 0 then begin    ; check ground state has all data - other states should even out
       if ms.oddeven[i] ne 'u' and ms.oddeven[i] ne 'g' and keyword_set(errorfile) then begin
           openw, lun, errorfile, /get_lun, /append
           printf, lun, ms.molid + 'bad odd even'
           close, lun
           free_lun, lun
       endif
       if ms.lambda[i] eq [0] and ms.plusmin[i] ne '+' and ms.plusmin[i] ne '-' and keyword_set(errorfile) then begin
           openw, lun, errorfile, /get_lun, /append
           printf, lun,  ms.molid + ': bad symm'
           close, lun
           free_lun, lun
       endif
       
    endif
    
    ; check for split states where 2S+1 is not 1
    if i ne n_elements(ms.name) - 1 then begin
    if  (ms.name[i+1] eq ms.name[i]) and (ms.smult[i] ne 1) and (ms.smult[i+1] ne 1) and keyword_set(errorfile) then begin
       openw, lun, errorfile, /get_lun, /append
       printf, lun,  ms.molid + ': duplicate state? : ' + ms.name[i] + ' --- ' + ms.name[i+1]
       close, lun
       free_lun, lun
    endif
    endif
    
    if keyword_set(showinfo) then print, 'Te = ', ms.Te[i], '   we = ', ms.we[i], '   Be = ', ms.Be[i]
    
    if ms.Te[i] eq 0. then grnd = 1
    
   ; if we don't know wxe, assume = 0.
   ; then there are infinite vibrational states
   ; but this shouldn't be a problem
   ; the sum converges fast anyway, so we just take a suitably large number
   
   if ms.wxe[i] ne 0. then begin
      vmax = fix(ms.we[i]/ms.wxe[i]/2. - 0.5)
   endif else begin
      vmax = 100
   endelse
   
   for v = 0, vmax do begin
   
      ; precompute some oft used terms 
      vph = v + 0.5
      
      ; compute Bv
      ; check for convergence
      t0 = ms.Be[i]
      t1 = ms.alfe[i]*vph
      t2 = ms.game[i]*vph*vph      
      if abs(t0) gt abs(t1*2) then begin
        BvST = t0 - t1
        if abs(t1) gt abs(t2*2) then begin      
           Bv = t0 - t1 + t2
        endif else begin
           Bv = t0 - t1
        endelse
      endif else begin
        Bv = t0
        BvST = t0
      endelse
 
      ; compute Dv     
      ; check for convergence       
      t0 = ms.De[i]
      t1 = ms.bete[i]*vph
      Dv = ms.De[i]
      if abs(t0) gt abs(t1*2) then begin
        Dv = t0 - t1
      endif else begin
        Dv = t0
      endelse
      
      if Bv le 0. then begin
         if ms.Te[i] eq 0. then begin
           if keyword_set(errorfile) then begin
             openw, lun, errorfile, /get_lun, /append
             printf, lun, ms.molid + ': ground state has Bv le 0'
             close, lun
             free_lun, lun
           endif
         endif
         goto, skip 
      endif
      
      ; compute G
      ; ensure series is asymptotic before including

      t0 = ms.we[i]*vph
      t1 = ms.wxe[i]*vph*vph
      t2 = ms.wye[i]*vph*vph*vph

      if abs(t0) gt abs(t1*2) then begin
        GST = t0 - t1
        if abs(t1) gt abs(t2*2) then begin      
           G = t0 - t1 + t2
        endif else begin
           G = t0 - t1
        endelse
      endif else begin
        G = t0
        GST = t0
      endelse
       
      ; calculate rotational part, explicitly
      
      Qrot = dblarr(n_elements(T))*0.
      QrotST = Qrot
      
      for k = 0, n_elements(T)-1 do begin
        Jmax = sqrt(1./(2*Bv*hckT[k]))    ; J where function is maximum
	    Jmax = Jmax * 10                  ; 10 times this maximum 
	    Jmax>=10                           ; very low T, Jmax -> 0
;    	Jmax<=1E2                         ; low or zero Bv
        if (Bv eq 0. and Dv eq 0.) then Jmax = 1 ; thus J=0 only and includes only first term, which should dominate 
        J = findgen(long(Jmax)) + ms.lambda[i]   ; J = Lambda, Lambda+1, ..... 
     
      ; calculate statistical weights for Lambda-doubling and nuclear structure
      ; see Irwin 1987
     
     ;print, ms.homonuc, ms.lambda[i], ms.plusmin[i], ms.oddeven[i]
     
     if ms.homonuc eq 0 then begin
         if ms.lambda[i] eq 0 then begin
            g_lh = 1
         endif else begin
            g_lh = 2
         endelse 
     endif else begin
         if ms.lambda[i] eq 0 then begin
            inuc = ms.I1
            gs = (inuc+1.)/(2.*inuc+1.)
            ga = inuc/(2.*inuc+1.)
            ieven = where(J mod 2 eq 0)
            g1 = J*0. + gs
            g1[ieven] = ga
            g2 = J*0. + ga
            g2[ieven] = gs
            g_lh = g1     ; we assume arbitrarily as default case
            if ((inuc mod 1) eq 0.) then begin   ; integral spin
               if ((ms.plusmin[i] eq '+' and ms.oddeven[i] eq 'g') or (ms.plusmin[i] eq '-' and ms.oddeven[i] eq 'u')) then g_lh = g2
            endif else begin
               if ((ms.plusmin[i] eq '-' and ms.oddeven[i] eq 'g') or (ms.plusmin[i] eq '+' and ms.oddeven[i] eq 'u')) then g_lh = g2
            endelse
         endif else begin
            g_lh = 1
         endelse 
     endelse
     
     
     ; total with fine structure (2S+1) 
     g_flh =  g_lh * ms.smult[i]
	 ;print, g_flh
	 
	 

        ; high order rotational terms can lead to divergence
	; ensure series is asymptotic before including
	
	F = Bv*J*(J+1)
	D = Dv*J*J*(J+1)*(J+1)
	ind = where(abs(F) gt abs(D*2), nind)
	if nind gt 0 then F[ind] = F[ind] - D[ind]
  
;        Qrot[k] = total((2*J+1)*exp(-hckT[k]*F)) / ms.sigma	
        qexp = -hckT[k]*(F+G+ms.Te[i]-Ezero)
        qexp2 = exp(qexp)
        ind3 = where(finite(qexp2) eq 0, nind3)
        if nind3 gt 0 then qexp2(ind3) = 0.d0
        
        Qrot[k] = total(g_flh*(2*J+1)*qexp2) 
;        Qrot[k] = total(g_flh*(2*J+1)*exp(-hckT[k]*(F+G+ms.Te[i]))) 
        
        
	FST = BvST*J*(J+1)
  	    
  	    qexpST = -hckT[k]*(FST+GST+ms.Te[i]-EzeroST)
        qexp2ST = exp(qexpST)
        ind3 = where(finite(qexp2ST) eq 0, nind3)
        if nind3 gt 0 then qexp2ST(ind3) = 0.d0
  	    
        QrotST[k] = total(g_flh*(2*J+1)*qexp2ST) 
        ;QrotST[k] = total(g_flh*(2*J+1)*exp(-hckT[k]*(FST+GST+ms.Te[i]))) 
      ;if ms.molid eq 'O2' and ms.Te[i] eq 0. then stop; print, g_flh  
       ; if ms.Te(i) eq 0 then stop
      endfor    
     
      Qtot = Qtot + Qrot ;* elvib_part
      if ms.Te[i] lt 40000. then Qtot4 = Qtot4 + Qrot
      if ms.Te[i] lt 40000. then QtotST = QtotST + QrotST
      
 
   endfor   
   skip:
   
endfor

;print, Qtot
;stop
Qi = Qtot ;* exp(hckT*Ezero)    ; Ezero is put above to avoid numerical problems at low T
Qi4 = Qtot4 ;* exp(hckT*Ezero)
QiST = QtotST ;* exp(hckT*EzeroST)

skipall:

if grnd eq 0 then begin
   print, 'No ground state! ns = ', ms.ns
   Qi = Qtot *0.
   Qi4 = Qtot4 *0.
   QiST = QtotST *0.
   if keyword_set(errorfile) then begin
   openw, lun, errorfile, /get_lun, /append
   printf, lun, ms.molid + ': no ground state  ns = ', ms.ns
   close, lun
   free_lun, lun
   endif
endif

   
end
