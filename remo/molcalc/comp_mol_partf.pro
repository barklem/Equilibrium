;
pro comp_mol_partf, T, molstruct, Qi

ms = molstruct ; for brevity, but retaining clarity in input line

h = 6.626076d-27   ; cgs
c = 2.997924d10
k = 1.38066d-16

hckT = h*c/k/T

Gzero = (0.5*ms.we[0] - 0.25*ms.wxe[0] + 0.125*ms.wye[0])   ;vibrational energy of lowest state in wavenumbers
fac1 = exp(hckT * Gzero)

Qtot = 0.d0
Qe = 0.d0

for i = 0, ms.ns-1 do begin

   ;if ms.Te[i] gt 200. then goto, skip  ; for stricter comparison
   
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
      
;      if Bv le 0. then goto, skip
      
      ; the electronic and vibrational parts in the summation
      ; putting fac1 inside here helps avoid machine limits
      
      if ms.lambda[i] eq 0 then begin
         gel = ms.smult[i]
      endif else begin
         gel = 2.*ms.smult[i]
      endelse 
      
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

;      elvib_part = gel * exp( -hckT * (G + ms.Te[i]) +hckT * Gzero)
      
;      ; calculate rotational part, in high T approximation
;      Qrot = 1.d0 / ms.sigma / hckT / Bv
      
      ; calculate rotational part, explicitly
      
      Qrot = dblarr(n_elements(T))*0.
      
      for k = 0, n_elements(T)-1 do begin
        Jmax = sqrt(1./(2*Bv*hckT[k]))    ; J where function is maximum
	Jmax = Jmax * 10                  ; 10 times this maximum 
	Jmax>=1                           ; very low T, Jmax -> 0
;	Jmax<=1E2                         ; low or zero Bv
        J = findgen(long(Jmax))           
	
	;print, Jmax

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
        Qrot[k] = total(gel*(2*J+1)*exp(-hckT[k]*(F+G+ms.Te[i]-Gzero))) $
	          / ms.sigma
      endfor    
     
      Qtot = Qtot + Qrot ;* elvib_part
 
      skip:
 
   endfor
   
endfor

Qi = Qtot
end
