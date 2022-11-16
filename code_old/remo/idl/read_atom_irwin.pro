pro read_atom_irwin, fname=fname, $
                     numco=numco, ncoef=ncoef, natom=natom, partfco=partfco,$
                     indx=indx, code=code, label=label, newlabel=newlabel, $
                     eion=eion

; reads atomic partition function fit coefficients from Irwin dataset
 
  
  if (n_elements(fname) eq 0) then $
     fname = '~/People/barklem/Equilibrium/remo/data/Irwin_partf/IRWIN_atoms.dat'
  
  maxatom = 1000L
  maxion  = 10L
  maxcoef = 10L

  indx    = lonarr(maxatom)
  code    = fltarr(maxatom)
  label   = strarr(maxatom)
  eion    = fltarr(maxatom)
  numco   = lonarr(maxatom)
  partfco = dblarr(maxcoef,maxatom)
  
  text     = ''
  newlabel = label

  
; open data file for reading
  openr, lun, /get_lun, fname

  
  i = 0L                        ; init counter
  while ~eof(lun) do begin
     readf,lun,text
     if (strmid(text,0,1) ne '#') then begin
        s  = strsplit(text,"'",/extract)
        s0 = strsplit(s[0]," ",/extract)
        s1 = strcompress(s[1],/remove_all)
        s2 = strsplit(s[2]," ",/extract)
        
        indx[i]	  = long(s0[0])
        code[i]	  = float(s0[1])
        label[i]  = s1
        eion[i]	  = float(s2[0])
        numco[i]  = long(s2[1])
        nco       = numco[i]
        partfco[0:nco-1,i] = double(s2[2:nco+1])
        
        len       = 0
        pos       = 0
        el        = ''        
        pos       = stregex(s1,'\++',len=len)
        if (pos ne -1) then begin
           el  = strmid(s1,0,pos)
           el  = strupcase(el)+'_'
           case len of
              1: el = el + 'II'
              2: el = el + 'III'
              3: el = el + 'IV'
              4: el = el + 'V'
              5: el = el + 'VI'
              6: el = el + 'VII'
           endcase
        endif else begin
           el = strupcase(s1)+'_I'
        endelse	
        newlabel[i] = el
        i++
     endif
  endwhile  

  natom   = i
  ncoef   = max(numco[0:natom-1])

  indx    = reform(indx[0:natom-1])
  code    = reform(code[0:natom-1])
  label   = reform(label[0:natom-1])
  eion    = reform(eion[0:natom-1])
  numco   = reform(numco[0:natom-1])
  partfco = reform(partfco[0:ncoef-1,0:natom-1],ncoef,natom)
  
  free_lun, lun

end
