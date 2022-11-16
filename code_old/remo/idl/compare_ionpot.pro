; script: compare ionisation potentials from CRC Handbook of Chemistry
; and Physics,  91st edition, and NIST
main_dir = '~/Dropbox/Shared/Equilibrium'
dir      = main_dir+'/remo/data/ionpot'
ionfile_crc   = 'ionpot_CRC_I-III.txt'
ionfile1_nist = 'ionpot1_nist.txt'
ionfile2_nist = 'ionpot2_nist.txt'
ionfile3_nist = 'ionpot3_nist.txt'

outdir  = main_dir+'/remo/output'
outfile = 'ion_table.txt'

; main program starts here
cd, current=current_dir

; cd to work directory
cd, dir

; read CRC data
readcol, ionfile_crc, atomnum_crc, elemsym_crc, ionpot_crc1, ionpot_crc2, ionpot_crc3, $
         format='I,A,F,F,F'

; read NIST data
read_ionpot_nist, ionfile1_nist, atomnum_arr1, elemsym_arr1, ionsym_arr1, ionpot_arr1, $
                  is_theor=is_theor1, is_interp=is_interp1
read_ionpot_nist, ionfile2_nist, atomnum_arr2, elemsym_arr2, ionsym_arr2, ionpot_arr2, $
                  is_theor=is_theor2, is_interp=is_interp2
read_ionpot_nist, ionfile3_nist, atomnum_arr3, elemsym_arr3, ionsym_arr3, ionpot_arr3, $
                  is_theor=is_theor3, is_interp=is_interp3

; print table with NIST results
openw,lunout,/get_lun,outdir+'/'+outfile
fmt_header ='(A1,X,A9,2X,A5,5X,A4,9X,A4,8X,A4)'
fmt_entries='(4X,I3,7X,A2,5X,F8.4,4X,F8.3,4X,F8.3)'
;printf, lunout, format=fmt_header, '*','Atom.Num.','Elem.','I','II','III'
printf, lunout, format=fmt_header, '*','Atom.Num.','Elem.','IE1 ','IE2 ','IE3 '
printf, lunout, format=fmt_header, '*','','','[eV]','[eV]','[eV]'
for iel=1L,92L do begin
   i = iel-1
   elemSym = elemSym_arr1[i]
   ionpot1 = ionpot_arr1[i]
   if iel eq 1 then begin
      ionpot2 = -1.0
      ionpot3 = -1.0
   endif else begin
      ionpot2 = ionpot_arr2[i-1]
      if iel eq 2 then begin
         ionpot3 = -1.0
      endif else begin
         ionpot3 = ionpot_arr3[i-2]
      endelse
   endelse
   
   printf, lunout, format=fmt_entries,iel, elemSym, ionpot1, ionpot2, ionpot3
endfor

free_lun,lunout
; end
cd, current_dir
end

