molec  = ['H2', 'H2+', 'H2O', 'OH', 'CH',  'CO', 'CN', 'C2', 'N2', 'O2', 'NO', 'NH']
nmolec = size(molec, /n_elem)
ncoef  = 5L
tabco  = fltarr(ncoef,nmolec)

for i=0,nmolec-1L do begin
 compare_lgkp, molec[i],co=co
 tabco[*,i]=co
endfor
end
