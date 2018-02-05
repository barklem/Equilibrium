;
pro comp_mol_partf, T, molstruct, Qi

; if keyword is set, reasonable estimates are made when we and Be are unknown

ms = molstruct ; for brevity, but retaining clarity in input line

h = 6.626076d-27   ; cgs
c = 2.997924d10
k = 1.38066d-16

hckT = h*c/k/T

; assume data are ordered highest to lowest
ilow = ms.ns-1
Gzero = (0.5*ms.we[ilow] - 0.25*ms.wxe[ilow] + 0.125*ms.wye[ilow])   ;vibrational energy of lowest state in wavenumbers
fac1 = exp(hckT * Gzero)

Qtot = 0.d0
Qe = 0.d0

; fill some data
ind = where (ms.we eq 0., nind, complement=ind2, ncomplement=nind2)
if nind ge 0 and nind2 ge 0 then ms.we[ind] = mean(ms.we[ind2]) 
ind = where (ms.Be eq 0., nind, complement=ind2, ncomplement=nind2)
if nind ge 0 and nind2 ge 0 then ms.Be[ind] = mean(ms.Be[ind2]) 

for i = 0, ms.ns-1 do begin
    
    print, 'Te = ', ms.Te[i], '   we = ', ms.we[i], '   Be = ', ms.Be[i]

    if ms.Te[i] lt 0. then goto, skip  ; no energy level
    
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
        if abs(t1) gt abs(t2*2) then begin      
           Bv = t0 - t1 + t2
        endif else begin
           Bv = t0 - t1
        endelse
      endif else begin
        Bv = t0
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
      
      if Bv le 0. then goto, skip
      
      ; the electronic and vibrational parts in the summation
      ; putting fac1 inside here helps avoid machine limits
      
      
      ; compute G
      ; ensure series is asymptotic before including

      t0 = ms.we[i]*vph
      t1 = ms.wxe[i]*vph*vph
      t2 = ms.wye[i]*vph*vph*vph

      if abs(t0) gt abs(t1*2) then begin
        if abs(t1) gt abs(t2*2) then begin      
           G = t0 - t1 + t2
        endif else begin
           G = t0 - t1
        endelse
      endif else begin
        G = t0
      endelse
       
      ; calculate rotational part, explicitly
      
      Qrot = dblarr(n_elements(T))*0.
      
      for k = 0, n_elements(T)-1 do begin
        Jmax = sqrt(1./(2*Bv*hckT[k]))    ; J where function is maximum
	Jmax = Jmax * 10                  ; 10 times this maximum 
	Jmax>=1                           ; very low T, Jmax -> 0
;	Jmax<=1E2                         ; low or zero Bv
        J = findgen(long(Jmax))  
     
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
            g_lh = g1
            if inuc eq inuc mod 1 then begin   ; integral spin
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
	
;	ind = where(F lt 0, nind)
;	if nind gt 0 then print, J[ind[0]], v, J[ind[0]]/Jmax 
	
;	print, J,F,	Bv*J*(J+1),Dv*J*J*(J+1)*(J+1)		  
;        Qrot[k] = total((2*J+1)*exp(-hckT[k]*F)) / ms.sigma		  
        Qrot[k] = total(g_flh*(2*J+1)*exp(-hckT[k]*(F+G+ms.Te[i]-Gzero))) 
        ;if k eq 0 then print, k, Qrot[k], F+G, ms.Te[i]-Gzero
      endfor    
     
      Qtot = Qtot + Qrot ;* elvib_part
      skip:
 
   endfor
   
endfor

Qi = Qtot
end
