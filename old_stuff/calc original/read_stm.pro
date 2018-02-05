pro read_stm, molid, A, B, D

; reads Sauval & Tatum data for molecules

 openr, lun, 'ST_molecules.dat', /get_lun

molid = strarr(500)
D = dblarr(500)
A = dblarr(5, 500)
B = dblarr(6, 500)

ns = 0
startp:

  s = ' '
  readf, lun, s, format='(a150)'
  if strmid(s,0,3) eq 'end' then goto, endp
  if strmid(s,0,1) eq '#' then goto, startp

  s = strsplit(s,' ', /extract)
  molid[ns] = s[0]
  D[ns] = double(s[1])
  B[0,ns] = double(s[2])
  B[1,ns] = double(s[3])
  B[2,ns] = double(s[4])
  B[3,ns] = double(s[5])
  B[4,ns] = double(s[6])
  B[5,ns] = double(s[7])
  A[0,ns] = double(s[8])
  A[1,ns] = double(s[9])
  A[2,ns] = double(s[10])
  A[3,ns] = double(s[11])
  A[4,ns] = double(s[12])

  ns = ns + 1
    
  goto, startp
  endp:
  
free_lun, lun

molid = molid[0:ns-1]
D = D[0:ns-1]
A = A[*,0:ns-1]
B = B[*,0:ns-1]

return
end
