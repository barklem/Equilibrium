pro convert_nist_table, fname, fnameout
  
  openr, lun, /get_lun, fname
  openw, lunout, /get_lun, fnameout

  s = ''

  printf, lunout,  '    J |  E                     |'
  printf, lunout,  ' -----|------------------------|'

  while ~eof(lun) do begin

     readf, lun, s
     has_data = ~( strmatch( s, '*Images*') )
     if (has_data) then begin
        res = strsplit( s, '<td>', /regex, /extract )
        s1  = res[3]
        s1a = strjoin(strsplit(s1,'</td>', /regex, /extract))
        s1b = strjoin(strsplit(s1a,'&nbsp;', /regex, /extract))
        if (strcompress(s1b,/remove_all) eq '') then s1b='---'
        s2  = res[4]
        s2a = strjoin(strsplit(s2,'</td>', /regex, /extract))
        s2b = strjoin(strsplit(s2a,'&nbsp;', /regex, /extract))
        s2c = strjoin(strsplit(s2b,'<b>', /regex, /extract))
        s2d = strjoin(strsplit(s2c,'</b>', /regex, /extract))
        s2e = strjoin(strsplit(s2d,'\?', /regex, /extract))
        printf, lunout, format='(X,A4,X,A,3X,A-21,A)', s1b, '|', s2e,'|'
        ;stop
     endif
     
  endwhile

  free_lun, lun, lunout

end
