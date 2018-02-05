dir     = '~/Desktop'
infile  = dir +'/'+ 'atomic_weights_NIST.txt'
outfile = dir +'/'+ 'standard_atomic_weights.txt'

dum = ''

openw, lun_out, outfile, /get_lun

openr, lun_in, infile, /get_lun

while ~eof(lun_in) do begin
      readf, lun_in, dum
      is_separator = stregex(dum,'_+',/bool) 
      
      if is_separator then begin
         if ~eof(lun_in) then begin 
            readf,  lun_in, dum
            printf, lun_out, dum
         endif
      endif
;stop
endwhile

free_lun, lun_in, lun_out

end
