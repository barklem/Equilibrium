pro reord_mol, ms, outstruct

; takes structure ms and reorders states in new strutcure msn


nstr = n_elements(ms.label)  ; biggest number of states is 46, make all structures dimension 50

label = strarr(nstr)
Te_str = replicate(' .', nstr)
w_str = strarr(nstr)
wx_str = strarr(nstr)
wy_str = strarr(nstr)
B_str = strarr(nstr)
alf_str = strarr(nstr)
gam_str = strarr(nstr)
D_str = strarr(nstr)
bet_str = strarr(nstr)
re_str = strarr(nstr)

Te_ref = strarr(nstr)
we_ref = strarr(nstr)
wxe_ref = strarr(nstr)
wye_ref = strarr(nstr)
Be_ref = strarr(nstr)
alfe_ref = strarr(nstr)
game_ref = strarr(nstr)
De_ref = strarr(nstr)
bete_ref = strarr(nstr)
re_ref = strarr(nstr)

nm = strarr(nstr)
lmbda = intarr(nstr)
spn = fltarr(nstr)
plsmn = strarr(nstr)
evenodd = strarr(nstr)

Te = double(clean_Te(ms.Te))   ; create a cleaned version, numbers only, unidentified are -1
Te_start = Te
ind = where(Te_start lt 0., nind)    ; move any states with no Te to the top
if nind gt 0 then Te[ind] = 1./(ind+1.)*1e15   ; arbitrary large number, trying to retain original ordering
ind = where(Te_start lt 0. and (ms.name eq '' and ms.label eq ''), nind)  ; move unlabelled states to bottom
if nind gt 0 then Te[ind] = -20.   ; arbitrary large number negative

ind2 = reverse(sort(Te))
label = ms.label[ind2]
nm = ms.name[ind2]
lmbda = ms.lambda[ind2]
spn = ms.spin[ind2]
plsmn = ms.plusmin[ind2]
evenodd = ms.oddeven[ind2]
Te_str = ms.Te[ind2]
w_str = ms.we[ind2]
wx_str = ms.wxe[ind2]
wy_str = ms.wye[ind2]
B_str = ms.Be[ind2]
alf_str = ms.alfe[ind2]
gam_str = ms.game[ind2]
D_str = ms.De[ind2]
bet_str = ms.bete[ind2]
re_str = ms.re[ind2]
Te_ref = ms.Te_ref[ind2]
we_ref = ms.we_ref[ind2]
wxe_ref = ms.wxe_ref[ind2]
wye_ref = ms.wye_ref[ind2]
Be_ref = ms.Be_ref[ind2]
alfe_ref = ms.alfe_ref[ind2]
game_ref = ms.game_ref[ind2]
De_ref = ms.De_ref[ind2]
bete_ref = ms.bete_ref[ind2]
re_ref = ms.re_ref[ind2]

;Te_new = double(clean_Te(Te_str))
;ind3 = where(Te gt 1e10, nind3)
;if nind3 gt 0 then Te_str[ind3] = replicate(' .', nind3)

outstruct = { $
    label  : label,       $
    name   : nm,          $
    lambda : lmbda,       $
    spin   : spn,         $
    plusmin: plsmn,       $
    oddeven: evenodd,       $
    Te    : Te_str,      $
    we   :    w_str,        $
    wxe  :    wx_str,       $
    wye  :    wy_str,       $
    Be   :    B_str,        $
    alfe :    alf_str,      $
    game :    gam_str,      $
    De   :    D_str,        $
    bete :    bet_str,      $
    re   :    re_str,         $
    Te_ref   :  Te_ref,      $
    we_ref   :  we_ref ,        $
    wxe_ref  : wxe_ref  ,       $
    wye_ref  : wye_ref  ,       $
    Be_ref   : Be_ref   ,        $
    alfe_ref : alfe_ref  ,      $
    game_ref : game_ref  ,      $
    De_ref   : De_ref  ,        $
    bete_ref : bete_ref  ,      $
    re_ref   : re_ref           $
} 
end