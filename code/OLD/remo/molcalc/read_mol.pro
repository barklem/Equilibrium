; reads in the molecular data into a structure
;
; Paul Barklem Feb 2007
;
pro read_mol, molfile, molstruct, debug=debug 

openr, lun, molfile, /get_lun

ns = 0
s = ' '
readf, lun, s, format='(a150)'; comment
readf, lun, s, format='(a150)'
s = strsplit(s,' ', /extract)
D = double(s[0])
sigma = fix(s[1])
mu = double(s[2])
ip = 0.d0
if n_elements(s) gt 3 then ip = double(s[3])
s=' '
readf, lun, s, format='(a150)' ; comment line

state = strarr(100)
Te = dblarr(100)
smult = intarr(100)
lambda = intarr(100)
Be = dblarr(100)
alfe = dblarr(100)
game = dblarr(100)
we = dblarr(100)
wxe = dblarr(100)
wye = dblarr(100)
De = dblarr(100)*0.d0
bete = dblarr(100)*0.d0
re = dblarr(100)*0.d0

ns = 0
debug1 = 0
if keyword_set(debug) then debug1 = 1

for i = 0, 99 do begin
  s = ' '
if debug1 eq 1 then print, '----'
  readf, lun, s, format='(a150)'
if debug1 eq 1 then print, s
  if strmid(s,0,3) eq 'end' then goto, breakout
  if strmid(s,0,1) eq '#' then goto, comment
  s = strsplit(s,' ', /extract)
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
  if n_elements(s) gt 10 then begin
  print, s[10]
  print, s[11]
  print, s[12]
  endif
endif
  state[ns] = s[0]
  smult[ns] = fix(s[1])
  lambda[ns] = fix(s[2])
  Te[ns] = double(s[3])
  we[ns] = double(s[4])
  wxe[ns] = double(s[5])
  wye[ns] = double(s[6])
  Be[ns] = double(s[7])
  alfe[ns] = double(s[8])
  game[ns] = double(s[9])
  if n_elements(s) gt 10 then begin
  De[ns] = double(s[10])
  bete[ns] = double(s[11])
  re[ns] = double(s[12])
  endif
  ns = ns + 1
  comment:
endfor

breakout:

free_lun, lun

state = state[0:ns-1]
smult = smult[0:ns-1]
lambda = lambda[0:ns-1]
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

molstruct = { $
    ns   :    ns,        $
    sigma:    sigma,     $
    D    :    D,         $
    ip   :    ip,        $
    mu   :    mu,        $
    state:    state,     $
    smult:    smult,     $
    lambda:   lambda,    $
    Te   :    Te,        $
    we   :    we,        $
    wxe  :    wxe,       $
    wye  :    wye,       $
    Be   :    Be,        $
    alfe :    alfe,      $
    game :    game,      $
    De   :    De,        $
    bete :    bete,      $
    re   :    re         $
  }


end
