; reads in the molecular data into a structure
;
; Paul Barklem Feb 2007
;              Jun-Sep 2010, new format
;

function clean_num, s
; cleans out any unwanted stuff around a number

snew = repstr(s, '(', '')
snew = repstr(snew, ')', '')
snew = repstr(snew, '[', '')
snew = repstr(snew, ']', '')
snew = repstr(snew, ' .', '0.0')   ; indicates no data
return, snew
end

function clean_Te, s
; cleans out any unwanted stuff around a number

snew = repstr(s, '(', '')
snew = repstr(snew, ')', '')
snew = repstr(snew, '[', '')
snew = repstr(snew, ']', '')
snew = repstr(snew, ' .', '-1.0')   ; indicates no data
return, snew
end

; main routine

pro read_mol, molfile, molstruct, debug=debug 

openr, lun, molfile, /get_lun

ns = 0
s = ' '
readf, lun, s, format='(a150)' ; comment
readf, lun, s, format='(a150)'
s = strsplit(s,' ', /extract)
molid = s[0]
I1 = float(s[1])
I2 = float(s[2])
homonuc = 0
if strpos(molid, '2') ge 0 then homonuc = 1
s = ' '
readf, lun, s, format='(a150)' ; comment
readf, lun, s, format='(a150)'
s = strsplit(s,' ', /extract)
D = double(s[0])
mu = double(s[1])
ip = 0.d0
if n_elements(s) gt 2 then ip = double(s[2])
s=' '
readf, lun, s, format='(a150)' ; comment line

state = strarr(100)
label = strarr(100)
lambda = intarr(100)
smult = intarr(100)
plusmin = strarr(100)
oddeven = strarr(100)
Te = dblarr(100)
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
  readf, lun, s, format='(a500)'
if debug1 eq 1 then print, s
  if strmid(s,0,3) eq 'end' then goto, breakout
  if strmid(s,0,1) eq '#' then goto, comment
  snew = strarr(16)
  reads, s, snew, format = '(a40, a7, a3, a3, a3, a3, 10a12)'
  s = snew
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
  state[ns] = s[0]
  label[ns] = s[1]
  lambda[ns] = fix(s[2])
  smult[ns] = fix(s[3])
  plusmin[ns] = strtrim(s[4])
  oddeven[ns] = strtrim(s[5])
  Te[ns] = double(clean_Te(s[6]))
  we[ns] = double(clean_num(s[7]))
  wxe[ns] = double(clean_num(s[8]))
  wye[ns] = double(clean_num(s[9]))
  Be[ns] = double(clean_num(s[10]))
  alfe[ns] = double(clean_num(s[11]))
  game[ns] = double(clean_num(s[12]))
  De[ns] = double(clean_num(s[13]))
  bete[ns] = double(clean_num(s[14]))
  re[ns] = double(clean_num(s[15]))
  ns = ns + 1
  comment:
endfor

breakout:

close, lun
free_lun, lun

state = state[0:ns-1]
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

molstruct = { $
    ns   :    ns,        $
    molid:    molid,     $
    homonuc:  homonuc,   $
    I1   :    I1,        $
    I2   :    I2,        $
    D    :    D,         $
    ip   :    ip,        $
    mu   :    mu,        $
    state:    state,     $
    label:    label,     $
    lambda:   lambda,    $
    smult:    smult,     $
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
    re   :    re         $
  }

end
