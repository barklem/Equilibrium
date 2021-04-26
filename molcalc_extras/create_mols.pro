pro create_mols

; creates final molecule data files from various sources
;
; Paul Barklem 10 Oct 2011

hhdir = 'HH_MolData'   ; directory with H&H data, no slash

; load the list of molecules, dissociation energies and component atoms

read_diss, molname, DE, components=molcomp
nmol = n_elements(molname)

; read structure containing data from Irikura (2007)

read_irikura, irikura

; read nuclear spins

read_spins, sname, spins

; loop through all molecules

lastWe = 0.
lastBe = 0.

for imol = 0, nmol-1 do begin

   ; extract nuclear spins
   atom1 = molcomp[imol, 0]
   atom2 = molcomp[imol, 1]
   atom1 = repstr(atom1, '+','')
   atom1 = repstr(atom1, '-','')
   atom2 = repstr(atom2, '+','')
   atom2 = repstr(atom2, '-','')
   ind1 = where(sname eq atom1)
   ind2 = where(sname eq atom2)
   I1 = spins[ind1]
   I2 = spins[ind2]
   ;print, atom1, I1
   ;print, atom2, I2

   molstr = 0.
   moldat = 0.
   ; read the NIST (i.e. Huber & Herzberg) data

   infile = hhdir + '/' + molname[imol] + '.html'
   outfile = 'HH_Parsed/' + molname[imol] + '.dat'    ; store for inspection
   parse_mol, infile, molstr, numstruct = moldat, file = outfile
  
   molstr_hh = molstr                                 ; store for inspection
   
   ; replace ground state data with that from Irikura if available
   
   ind = where(molname[imol] eq irikura.name, nind)
   if nind gt 0 then begin                            ; molecule in list
       ind2 = where(moldat.Te eq 0.0, nind2)
       if nind2 gt 0 then begin
           molstr.we[ind2] = strcompress(string(irikura.we[ind]), /remove_all)
           molstr.wxe[ind2] = strcompress(string(irikura.wxe[ind]), /remove_all)
           molstr.wye[ind2] = strcompress(string(irikura.wye[ind]), /remove_all)
           molstr.Be[ind2] = strcompress(string(irikura.Be[ind]), /remove_all)
           molstr.alfe[ind2] = strcompress(string(irikura.alfe[ind]), /remove_all)
           molstr.we_ref[ind2] = replicate('IR', nind2)    ; update references
           molstr.wxe_ref[ind2] = replicate('IR', nind2)
           molstr.wye_ref[ind2] = replicate('IR', nind2)
           molstr.Be_ref[ind2] = replicate('IR', nind2)
           molstr.alfe_ref[ind2] = replicate('IR', nind2)
           print, 'Replaced ', molname[imol], ' ground state data with data from Irikura (2007)'
       endif
   endif
   

   
   ; look for additional replacement data
   ; how to do this?
   
   infile = 'Replacements/' + molname[imol] + '.rep'
   q = file_test(infile)
   if q eq 1 then begin
      read_rep, infile, repstr
      if repstr.ns ge 1 then begin
      for istate = 0, repstr.ns - 1 do begin
          ind = where (((repstr.name[istate] eq molstr.name)           $     ; look for states with same name and spectroscopic term
                        and (repstr.spin[istate] eq molstr.spin)      $
                        and (repstr.lambda[istate] eq molstr.lambda)) or $
                        ((repstr.Te[istate] eq molstr.Te))           $        ; look for states with same Te - allows names to be changed
                        , nind)                                              ; and thus data to be changed in two steps
          if nind gt 0 then begin   ; replace data
               
             inew = ind[0] 
             print, 'Some data replaced in ', molname[imol]
             
             
          endif else begin          ; add data to structure
             ind2 = where(molstr.name eq '' and moldat.Te le 0., nind2)
             inew = ind2[0]
             print, 'Some data added in ', molname[imol]
             if strcompress(string(repstr.name[istate]), /remove_all) ne '.' then molstr.label[inew] = strcompress(string(repstr.name[istate]), /remove_all)
             
             ; need blank to fill
             
             molstr.name[inew] = '.'
             molstr.lambda[inew] = -1
             molstr.spin[inew] =  -1
             molstr.plusmin[inew] = '.'
             molstr.oddeven[inew] = '.'
             molstr.Te[inew] = '.'
             molstr.we[inew] = '.'
             molstr.wxe[inew] = '.'
             molstr.wye[inew] = '.'
             molstr.Be[inew] = '.'
             molstr.alfe[inew] = '.'
             molstr.game[inew] = '.'
             molstr.De[inew] = '.'
             molstr.bete[inew] = '.'
             molstr.re[inew] = '.'
             molstr.Te_ref[inew] = '.'
             molstr.we_ref[inew] = '.'
             molstr.wxe_ref[inew] = '.'
             molstr.wye_ref[inew] = '.'
             molstr.Be_ref[inew] = '.'
             molstr.alfe_ref[inew] = '.'
             molstr.game_ref[inew] = '.'
             molstr.De_ref[inew] = '.'
             molstr.bete_ref[inew] = '.'
             molstr.re_ref[inew] = '.'
          endelse
             if repstr.chlab[istate] gt 0 then molstr.label[inew] = repstr.label[istate]
             if strcompress(string(repstr.name[istate]), /remove_all) ne '.' then molstr.name[inew] = strcompress(string(repstr.name[istate]), /remove_all)
             if repstr.lambda[istate] ge 0 then molstr.lambda[inew] = repstr.lambda[istate]
             if repstr.spin[istate] ge 0 then molstr.spin[inew] = repstr.spin[istate]
             if strcompress(string(repstr.plusmin[istate]), /remove_all) ne '.' then molstr.plusmin[inew] = strcompress(string(repstr.plusmin[istate]), /remove_all)
             if strcompress(string(repstr.oddeven[istate]), /remove_all) ne '.' then molstr.oddeven[inew] = strcompress(string(repstr.oddeven[istate]), /remove_all)
             if double(clean_Te(repstr.Te[istate])) ne -1.0 then molstr.Te[inew] = strcompress(repstr.Te[istate], /remove_all)
             if double(clean_num(repstr.we[istate])) ne 0. then molstr.we[inew] = strcompress(repstr.we[istate], /remove_all)
             if double(clean_num(repstr.wxe[istate])) ne 0. then molstr.wxe[inew] = strcompress(repstr.wxe[istate], /remove_all)
             if double(clean_num(repstr.wye[istate])) ne 0. then molstr.wye[inew] = strcompress(repstr.wye[istate], /remove_all)
             if double(clean_num(repstr.Be[istate])) ne 0. then molstr.Be[inew] = strcompress(repstr.Be[istate], /remove_all)
             if double(clean_num(repstr.alfe[istate])) ne 0. then molstr.alfe[inew] = strcompress(repstr.alfe[istate], /remove_all)
             if double(clean_num(repstr.game[istate])) ne 0. then molstr.game[inew] = strcompress(repstr.game[istate], /remove_all)
             if double(clean_num(repstr.De[istate])) ne 0. then molstr.De[inew] = strcompress(repstr.De[istate], /remove_all)
             if double(clean_num(repstr.bete[istate])) ne 0. then molstr.bete[inew] = strcompress(repstr.bete[istate], /remove_all)
             ; molstr.re[inew] = string(repstr.re[istate])   not used
             if double(clean_Te(repstr.Te[istate])) ne -1.0 then molstr.Te_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.we[istate])) ne 0. then molstr.we_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.wxe[istate])) ne 0. then molstr.wxe_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.wye[istate])) ne 0. then molstr.wye_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.Be[istate])) ne 0. then molstr.Be_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.alfe[istate])) ne 0. then molstr.alfe_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.game[istate])) ne 0. then molstr.game_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.De[istate])) ne 0. then molstr.De_ref[inew] = repstr.ref[istate]
             if double(clean_num(repstr.bete[istate])) ne 0. then molstr.bete_ref[inew] = repstr.ref[istate]
    endfor
    endif
   endif
   ;if molname[imol] eq 'CH+'then stop
   
   ; reorder if replacements have been made
  
   if q eq 1 then begin 
     reord_mol, molstr, molstr_new
     molstr = molstr_new
   endif
   
   
   ; check for cases where we for the ground state is not known
   ; and make an estimate based on nearby molecules
   
   tTe = double(clean_Te(molstr.Te))
   tWe = double(clean_num(molstr.We))
   ind = where(tTe eq 0., nind)
   ind2 = where(tWe gt 0., nind2)   ; if we exists for other states, this will be preferred
   if nind gt 0 then begin
      for ii = 0, nind-1 do begin
          if tWe[ind[ii]] le 0. and nind2 le 0 then begin 
             molstr.We[ind[ii]] = string(round(lastWe*10.)/10., '(f10.1)')   ; assign and round to 1 decimals
             molstr.We_ref[ind[ii]] = 'EST'
          endif
        ;  if tBe[ind[ii]] le 0. and nind2 gt 0 then begin 
        ;     molstr.Be[ind[ii]] = string(round(mean(tBe[ind2])*1000.)/1000.)   ; assign and round 
        ;     molstr.Be_ref[ind[ii]] = 'EST'
        ;     print, '********', molname[imol]
        ;endif
      endfor
   endif
   
   ; check for cases where Be for the ground state is not known
   ; and make an estimate based on nearby molecules
   
   tTe = double(clean_Te(molstr.Te))
   tBe = double(clean_num(molstr.Be))
   ind = where(tTe eq 0., nind)
   ind2 = where(tBe gt 0., nind2)   ; if Be exists for other states, this will be preferred
   if nind gt 0 then begin
      for ii = 0, nind-1 do begin
          if tBe[ind[ii]] le 0. and nind2 le 0 then begin 
             molstr.Be[ind[ii]] = string(round(lastBe*1000.)/1000., '(f10.3)')   ; assign and round to one decimal
             molstr.Be_ref[ind[ii]] = 'EST'
          endif
        ;  if tBe[ind[ii]] le 0. and nind2 gt 0 then begin 
        ;     molstr.Be[ind[ii]] = string(round(mean(tBe[ind2])*1000.)/1000.)   ; assign and round 
        ;     molstr.Be_ref[ind[ii]] = 'EST'
        ;     print, '********', molname[imol]
        ;endif
      endfor
   endif
   
   tWe = double(clean_num(molstr.We))
   tBe = double(clean_num(molstr.Be))
   
   ; fill some more data
   ind = where (tWe eq 0. and tTe ge 0., nind)
   ind2 = where(tWe gt 0., nind2)
   if nind gt 0 and nind2 gt 0 then begin
       molstr.we[ind] = string(round(mean(tWe[ind2])*10.)/10., '(f10.1)')   ; assign and round to one decimal
       molstr.we_ref[ind] = 'EST'  
   endif
   ind = where (tBe eq 0. and tTe ge 0., nind)
   ind2 = where(tBe gt 0., nind2)
   if nind gt 0 and nind2 gt 0 then begin
       molstr.Be[ind] = string(round(mean(tBe[ind2]) *1000.)/1000., '(f10.3)')   ; assign and round to one decimal
       molstr.Be_ref[ind] = 'EST'  
   endif
          
   ;if molname[imol] eq 'CBr' then stop

   
   ; make a structure with numbers, in format needed by codes
   
    ind = where(double(clean_Te(molstr.Te)) ge 0.d0 or molstr.label ne '', ns)
    homonuc = 0
    if strpos(molname[imol], '2') ge 0 then homonuc = 1
    
    moldat = { $
    ns     : ns,   $
    molid  : molname[imol], $
    I1 : I1[0], $
    I2 : I2[0], $
    D : DE[imol], $
    homonuc: homonuc, $
    label  : molstr.label,       $
    name   : molstr.name,          $
    lambda : molstr.lambda,       $
    smult  : molstr.spin,         $
    plusmin: molstr.plusmin,       $
    oddeven: molstr.oddeven,       $
    Te   :    double(clean_Te(molstr.Te)),      $
    we   :    double(clean_num(molstr.we)),        $
    wxe  :    double(clean_num(molstr.wxe)),       $
    wye  :    double(clean_num(molstr.wye)),       $
    Be   :    double(clean_num(molstr.Be)),        $
    alfe :    double(clean_num(molstr.alfe)),      $
    game :    double(clean_num(molstr.game)),      $
    De   :    double(clean_num(molstr.De)),        $
    bete :    double(clean_num(molstr.bete)),      $
    re   :    double(clean_num(molstr.re)),         $
    Te_ref   :  molstr.Te_ref,      $
    we_ref   :  molstr.we_ref ,        $
    wxe_ref  : molstr.wxe_ref  ,       $
    wye_ref  : molstr.wye_ref  ,       $
    Be_ref   : molstr.Be_ref   ,        $
    alfe_ref : molstr.alfe_ref  ,      $
    game_ref : molstr.game_ref  ,      $
    De_ref   : molstr.De_ref  ,        $
    bete_ref : molstr.bete_ref  ,      $
    re_ref   : molstr.re_ref           $
    }
    
    moldat_hh = { $
    ns     : ns,   $
    molid  : molname[imol], $
    I1 : I1[0], $
    I2 : I2[0], $
    D : DE[imol], $
    homonuc: homonuc, $
    label  : molstr_hh.label,       $
    name   : molstr_hh.name,          $
    lambda : molstr_hh.lambda,       $
    smult  : molstr_hh.spin,         $
    plusmin: molstr_hh.plusmin,       $
    oddeven: molstr_hh.oddeven,       $
    Te   :    double(clean_Te(molstr_hh.Te)),      $
    we   :    double(clean_num(molstr_hh.we)),        $
    wxe  :    double(clean_num(molstr_hh.wxe)),       $
    wye  :    double(clean_num(molstr_hh.wye)),       $
    Be   :    double(clean_num(molstr_hh.Be)),        $
    alfe :    double(clean_num(molstr_hh.alfe)),      $
    game :    double(clean_num(molstr_hh.game)),      $
    De   :    double(clean_num(molstr_hh.De)),        $
    bete :    double(clean_num(molstr_hh.bete)),      $
    re   :    double(clean_num(molstr_hh.re)),         $
    Te_ref   :  molstr_hh.Te_ref,      $
    we_ref   :  molstr_hh.we_ref ,        $
    wxe_ref  : molstr_hh.wxe_ref  ,       $
    wye_ref  : molstr_hh.wye_ref  ,       $
    Be_ref   : molstr_hh.Be_ref   ,        $
    alfe_ref : molstr_hh.alfe_ref  ,      $
    game_ref : molstr_hh.game_ref  ,      $
    De_ref   : molstr_hh.De_ref  ,        $
    bete_ref : molstr_hh.bete_ref  ,      $
    re_ref   : molstr_hh.re_ref           $
    }


   print, 'CHANGING HE2 ZERO POINT!!!!!!'
   if molname[imol] eq 'He2' then begin
    moldat.Te = moldat.Te - 144048.
    moldat_hh.Te = moldat_hh.Te - 144048.
    molstr.Te = string(moldat.Te)
    molstr_hh.Te = string(moldat_hh.Te)
   endif  
   
    ; save all data

    outfile1 = 'Final_Mols_idl/'+molname[imol]+'.idl'
    save, file = outfile1, moldat

    moldat_temp = moldat
    moldat = moldat_hh
    outfile1b = 'Final_Mols_idl/'+molname[imol]+'_hh.idl'
    save, file = outfile1b, moldat
    moldat = moldat_temp
    
    outfile2 = 'Final_Mols_ascii/'+molname[imol]+'.txt'
    D = DE[imol]
    openw, 1, outfile2
    printf, 1, 'D0', 'I1', 'I2', format = '(3a10)'
    printf, 1, D, I1[0], I2[0], format = '(f10.7, 2f10.1)'
    printf, 1, 'label', 'name', 'Lambda', '2S+1', '+/-', 'u/g', $
                'Te', 'we', 'wxe', 'wye', 'Be', 'alfe','game','De', 'bete', 're', 'references', $
                format = '(a40,a7,4a7,10a15,3x,a10)'
    printf, 1, ' ', ' ', ' ', ' ', ' ', ' ', $
                '[cm-1]', '[cm-1]', '[cm-1]', '[cm-1]', '[cm-1]', '[cm-1]','[cm-1]','[cm-1]', '[cm-1]', '[AA]', ' ', $
                format = '(a40,a7,4a7,10a15,3x,a10)'
    for j = 0, n_elements(molstr.label) - 1 do begin
       
       lambda_temp = string(molstr.lambda[j], '(i1)')
       if molstr.lambda[j] lt 0 then lambda_temp = '.' 
       spin_temp = string(molstr.spin[j], '(i2)')
       if molstr.spin[j] lt 0 then spin_temp = '.'
       
       if (molstr.label[j] ne '' or moldat.Te[j] ge 0.d0) then $
       printf, 1, molstr.label[j], molstr.name[j], lambda_temp, spin_temp, molstr.plusmin[j], $
       molstr.oddeven[j], molstr.Te[j], molstr.we[j], molstr.wxe[j], molstr.wye[j], molstr.Be[j], $
       molstr.alfe[j], molstr.game[j], molstr.De[j], molstr.bete[j], molstr.re[j], $
       molstr.Te_ref[j], molstr.we_ref[j], molstr.wxe_ref[j], molstr.wye_ref[j], $
       molstr.Be_ref[j], molstr.alfe_ref[j], molstr.game_ref[j], molstr.De_ref[j], $
       molstr.bete_ref[j], molstr.re_ref[j], format = '(a40, 5a7, 10a15, 10a7)'
    endfor
    printf, 1, 'END'
    close, 1

    
    ind2 = where(moldat.Te eq 0., nind2)
    if nind2 gt 0 then begin
       newlastBe = mean(moldat.Be[ind2])
       if newlastBe gt 0. then lastBe = newlastBe
    endif
    if nind2 gt 0 then begin
       newlastWe = mean(moldat.We[ind2])
       if newlastWe gt 0. then lastWe = newlastWe
    endif
    ;print, 'last', lastBe
    ;stop
    
   
;if molname[imol] eq 'FeH' then stop
endfor   
stop
end







