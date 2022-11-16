; script: extract reference codes from NIST_levels_ref file
; substitute cite-keys for codes from NIST_codes_key file

dir   = '~/Dropbox/Shared/Equilibrium/remo/references'
fname = 'NIST_levels_ref.txt'
keysfile = 'NIST_codes_keys.txt'
fout  = 'NIST_levels_ref_list.txt'

dumStr = ''

; current directory
cd, current=current_dir

; read codes and keys
cd, dir
readcol,keysfile,format='A,A',codeStrings,keyStrings

; open input file for reading
cd, dir
openr, lun, /get_lun, fname
;open output file for writing
openw, lunout, /get_lun, fout

; read header string
readf,lun,dumStr

; read lines, process them and rewrite them to output file
printf, lunout,'* Elem', 'Ion', 'References (citation keys)',format='(A,2X,A,4X,A)'
while ~EOF(lun) do begin
   readf,lun,dumStr
   dumStrArr = strsplit(dumStr,/extract)
   ndumStrArr = n_elements(dumStrArr)
   elemSym=dumStrArr[0]
   case dumStrArr[1] of
      '1': ionSym='I'
      '2': ionSym='II'
      '3': ionSym='III'
   endcase
   nrefs = ndumStrArr-2
   keyStrArr = strarr(nrefs)
   for iref=0,nrefs-1 do begin
      i = iref+2
      codeStr = dumStrArr[i]
      ww=where(codeStrings eq codeStr)
      keyStrArr[iref] = keyStrings[ww[0]]
   endfor
   printf,lunout,elemSym,ionSym,keyStrArr,format='(3X,A2,3X,A3,3X,7(X,A10))'   

endwhile

; close input and output files
free_lun, lun, lunout

; return to original directory
cd, current_dir
end
