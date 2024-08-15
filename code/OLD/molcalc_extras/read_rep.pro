; reads replacement data in simple form

pro read_rep, molfile, molstruct, debug=debug 

openr, lun, molfile, /get_lun

;s=' '
;readf, lun, s, format='(a150)' ; reference line
;reference = strcompress(s)
s=' '
readf, lun, s, format='(a150)' ; comment line

name = strarr(100)
label = strarr(100)   ; only changed by special command
lambda = intarr(100)
smult = intarr(100)
plusmin = strarr(100)
oddeven = strarr(100)
Te = strarr(100)
Be = strarr(100)
alfe = strarr(100)
game = strarr(100)
we = strarr(100)
wxe = strarr(100)
wye = strarr(100)
De = strarr(100)
bete = strarr(100)
re = strarr(100)
ref = strarr(100)
chlab = intarr(100)

ns = 0
debug1 = 0
if keyword_set(debug) then debug1 = 1

for i = 0, 99 do begin
  s = ' '
if debug1 eq 1 then print, '----'
  readf, lun, s, format='(a500)'
if debug1 eq 1 then print, s
  if strmid(s,0,3) eq 'end' then goto, breakout
  if strmid(s,0,1) eq '#' then goto, comment
  if strmid(s,0,12) eq 'change label' then begin
      label[ns-1] = strtrim(strmid(s,12,20))
      chlab[ns-1] = 1
      goto, comment
  endif
  s = strsplit(s, /extract)
  
if debug1 eq 1 then begin
  print, s[0]
  print, s[1]
  print, s[2]
  print, s[3]
  print, s[4]
  print, s[5]
  print, s[6]
  print, s[7]
  print, s[8]
  print, s[9]  
  print, s[10]
  print, s[11]
  print, s[12]
  print, s[13]
  print, s[14]
  print, s[15]
endif
;print, s
  name[ns] = s[0]
  lambda[ns] = fix(s[2])
  smult[ns] = fix(s[1])
  plusmin[ns] = strtrim(s[3])
  oddeven[ns] = strtrim(s[4])
  Te[ns] = strtrim(s[5])
  we[ns] = strtrim(s[6])
  wxe[ns] = strtrim(s[7])
  wye[ns] = strtrim(s[8])
  Be[ns] = strtrim(s[9])
  alfe[ns] = strtrim(s[10])
  game[ns] = strtrim(s[11])
  De[ns] = strtrim(s[12])
  bete[ns] = strtrim(s[13])
  re[ns] = strtrim(s[14])
  ref[ns] = strcompress(s[15])
  ns = ns + 1
  comment:
endfor

breakout:

close, lun
free_lun, lun

if ns ge 1 then begin

name = name[0:ns-1]
label = label[0:ns-1]
lambda = lambda[0:ns-1]
smult = smult[0:ns-1]
plusmin = plusmin[0:ns-1]
oddeven = oddeven[0:ns-1]
Te = Te[0:ns-1]
Be = Be[0:ns-1]
alfe = alfe[0:ns-1]
game = game[0:ns-1]
we = we[0:ns-1]
wxe = wxe[0:ns-1]
wye = wye[0:ns-1]
De = De[0:ns-1]
bete = bete[0:ns-1]
re = re[0:ns-1]
ref = ref[0:ns-1]

endif

molstruct = { $
    ns   :    ns,        $
    label:    label, $
    chlab:    chlab, $
    name:    name,     $
    lambda:   lambda,    $
    spin  :    smult,     $
    plusmin:  plusmin,   $
    oddeven:  oddeven,   $
    Te   :    Te,        $
    we   :    we,        $
    wxe  :    wxe,       $
    wye  :    wye,       $
    Be   :    Be,        $
    alfe :    alfe,      $
    game :    game,      $
    De   :    De,        $
    bete :    bete,      $
    re   :    re,        $
    ref  :    ref        $
  }

end
