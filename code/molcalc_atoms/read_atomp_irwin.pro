pro read_atomp_irwin, ind,code,label,eion,nco,co,label2,infile=infile,mlines=mlines,mco=mco

  if (n_elements(infile) eq 0) then infile = '/Users/remo/eqw/IRWIN_atoms.dat'
  if (n_elements(mlines) eq 0) then mlines = 300L
  if (n_elements(mco) eq 0) then mco = 6L
  

  text = ''

  ind	= lonarr(mlines)
  code  = fltarr(mlines)
  label = strarr(mlines)
  eion	= fltarr(mlines)
  nco	= lonarr(mlines)
  co	= dblarr(mco,mlines)

  label2 = strarr(mlines)

  openr, lun, infile, /get_lun


  i = 0L
  while ~eof(lun) do begin
     readf,lun,text
     if (strmid(text,0,1) ne '#') then begin
        s  = strsplit(text,"'",/extract)
        s0 = strsplit(s[0]," ",/extract)
        s1 = strcompress(s[1],/remove_all)
        s2 = strsplit(s[2]," ",/extract)
        
        ind[i]	 = long(s0[0])
        code[i]	 = float(s0[1])
        label[i] = s1
        eion[i]	 = float(s2[0])
        nco	 = long(s2[1])
        co[0:nco-1,i] = double(s2[2:nco+1])
        
        len  = 0
        pos  = 0
        el   = ''
        
        pos  = stregex(s1,'\++',len=len)
        if (pos ne -1) then begin
           el  = strmid(s1,0,pos)
           el  = strupcase(el)+'_' ;+'_I'
           case len of
              1: el=el+'II'
              2: el=el+'III'
              3: el=el+'IV'
              4: el=el+'V'
              5: el=el+'VI'
              6: el=el+'VII'
           endcase
        endif else begin
           el = strupcase(s1)+'_I'
        endelse	
        label2[i] = el
        i++
     endif
  endwhile  

  print, i
  free_lun, lun

end
