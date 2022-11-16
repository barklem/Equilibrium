pro read_sta, atomid, C, maxE

; reads Sauval & Tatum data for atoms

 openr, lun, 'ST_atom.dat', /get_lun

atomid = strarr(200)
maxE = dblarr(200)
C = dblarr(5, 200)

ns = 0
startp:

  s = ' '
  readf, lun, s, format='(a100)'
  if strmid(s,0,3) eq 'end' then goto, endp

  s = strsplit(s,' ', /extract)
  atomid[ns] = s[0]
  maxE[ns] = double(s[1])
  C[0,ns] = double(s[2])
  C[1,ns] = double(s[3])
  C[2,ns] = double(s[4])
  C[3,ns] = double(s[5])
  C[4,ns] = double(s[6])

  ns = ns + 1
    
  goto, startp
  endp:
  
free_lun, lun

atomid = atomid[0:ns-1]
maxE = maxE[0:ns-1]
C = C[*,0:ns-1]

return
end
