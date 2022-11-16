function extract_col, s, j
; extract column j from string s
sout = ' '
if j gt strcnt(s, '<td') then begin
   print, 'EXTRACT_COL WARNING: asks for column which does not exist'
   return, sout
endif
i = 0
cnt = 0
while (cnt ne j) do begin 
   i = strpos(s, '<td', i) 
   if (i ne -1) then begin  
     cnt = cnt + 1 
     i = i + 1  
   endif 
endwhile
ps = strpos(s, '>', i)+1
pe = strpos(s, '</td', ps)
sout = strmid(s, ps, pe-ps)
return, sout
end


function remove_links, s
snew = s
ps = 0
while ps ge 0 do begin
   ps = strpos(snew, '<a ')
   if ps ge 0 then begin
      pe = strpos(snew, '</a', ps) + 4
      sr = strmid(snew, ps, pe-ps)
      snew = repstr(snew, sr, '')
   endif
endwhile
return, snew
end

function clean_number, s
; cleans out any unwanted stuff around a number

; NOTE, order can be important!

iev = strpos(s, 'eV')

snew = repstr(s, '<sub>', '')
snew = repstr(snew, '</sub>', '')
snew = repstr(snew, '&nbsp;', '.')     ; blank
snew = repstr(snew, 'H<sup>Q</sup>', '')
snew = repstr(snew, 'H<sup>R</sup>', '')
snew = remove_links(snew)              ; must be above below, or else a's are removed
snew = repstr(snew, '(D2)', '')
snew = repstr(snew, 'B1 = ', '')
snew = repstr(snew, 'D1 = ', '')
snew = repstr(snew, 'r1 = ', '')
snew = repstr(snew, 'B1=', '')
snew = repstr(snew, 'D1=', '')
snew = repstr(snew, 'r1=', '')
snew = repstr(snew, 'B2 = ', '')
snew = repstr(snew, 'D2 = ', '')
snew = repstr(snew, 'r2 = ', '')
snew = repstr(snew, 'B1', '')
snew = repstr(snew, 'D1', '')
snew = repstr(snew, 'r1', '')
snew = repstr(snew, 'B2', '')
snew = repstr(snew, 'D2', '')
snew = repstr(snew, 'r2', '')
snew = repstr(snew, '~', '')
snew = repstr(snew, 'H', '')
snew = repstr(snew, 'Z', '')
snew = repstr(snew, 'R', '')
snew = repstr(snew, 'V', '')
snew = repstr(snew, 'i', '')
snew = repstr(snew, 'e', '')
snew = repstr(snew, '&gt;', '')
snew = repstr(snew, '&lt;', '')
snew = repstr(snew, '&#916;G(3/2) =', '')  ; Delta G(3/2) in SeN
snew = repstr(snew, '&#8804;', '')  ; leq
snew = repstr(snew, '&#8805;', '')  ; geq
;snew = repstr(snew, '*', '')
snew = strcompress(snew, /remove_all)
if iev ge 0 then snew = snew + 'eV'
if snew eq '' then snew = '.'

snew = repstr(snew, 'O.OO664', '0.00664') ; fixes a specific error in BeS state

return, snew
end

pro get_symmetry, label, name, Lambda, spin, plusmin, oddeven
; takes the processed label and extracts the symmetries 
; Lambda, 2S+1, +/-, g/u and the state name (e.g. X, a etc)
; Lambda and 2S+1 are returned -1 if not found
; - since +/- is only relevant for Sigma states
;   '-' = -, '+' = +, '.' = both in case of Lambda > 0, or undefined
; - since g/u symmetry is only relevant for homonuclear molecules
;   'g' = g, 'u' = u, '.' = not applicable/undefined

name = ' '
i = strpos(label, '^')
if i ge 0 then begin
   name = strmid(label, 0, i)
   name = repstr(name, '(','')
endif else begin
   i = strpos(label, '(')
   if i ge 0 then begin
      name = strmid(label, 0, i)
   endif else begin
      name = label
   endelse
endelse
i = strpos(name, 'Ome')
if i ge 0 then name = strmid(name, 0, i)
name = strtrim(name)
name = repstr(name, ')', '')
name = strcompress(name, /remove_all)

lampos = -1
Lambda = -1
i = strpos(label, 'Sigma')
if i ge 0 then begin
  Lambda = 0
  lampos = i
endif
i = strpos(label, 'Pi')
if i ge 0 then begin
  Lambda = 1
  lampos = i
endif
i = strpos(label, 'Delta')
if i ge 0 then begin
  Lambda = 2
  lampos = i
endif
i = strpos(label, 'Phi')
if i ge 0 then begin
  Lambda = 3
  lampos = i
endif
spin = -1
if lampos gt 0 then begin
  spinpos = strpos(label, '^', lampos, /reverse_search) + 1
  if spinpos gt 0 then begin
    spinstr = strmid(label, spinpos, 2)
    ;print, spinstr
    spin = float(spinstr)    
  endif  
endif

plusmin = '.'
if Lambda eq 0 then begin
  i = strpos(label, '^+', lampos)
  if i ge 0 then plusmin = '+'
  i = strpos(label, '^-', lampos)
  if i ge 0 then plusmin = '-'
  i = strpos(label, '^(+)', lampos)
  if i ge 0 then plusmin = '+'
  i = strpos(label, '^(-)', lampos)
  if i ge 0 then plusmin = '-'
endif

oddeven = '.'
i = strpos(label, '_g', lampos)
if i ge 0 then oddeven = 'g'
i = strpos(label, '_u', lampos)
if i ge 0 then oddeven = 'u'

end





function strarrcnt, s, ss
; return number of occurnences of ss in each element of array s
n = n_elements(s)
cnt = intarr(n)*0
for i = 0, n-1 do cnt[i] = strcnt(s[i], ss)
return, cnt
end

 
function strcnt, s, ss
; return number of occurences of ss in s
; s is a single string
i = 0 
cnt = 0 
while (i ne -1) do begin 
   i = strpos(s, ss, i) 
   if (i ne -1) then begin  
     cnt = cnt + 1 
     i = i + 1  
   endif 
endwhile
return, cnt
end


pro parse_mol, infile, outstruct, numstruct=numstruct, file=file, ion=ion

; parses html output from NIST
; obtained selecting:
;   - Gas phase ion energetics data
;   - Constants of diatomic molecules
;
;
; outstruct is strings with all info
; numstruct is numbers
; Paul Barklem June 2011


nstr = 50  ; biggest number of states is 46, make all structures dimension 50

label = strarr(nstr)
Te_str = replicate('-1.0', nstr)
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



q = file_test(infile)
if q ne 1 then goto, skip ; create empty structure

; read html file into string array
nl = file_line(infile)
s = strarr(nl)
openr, lun, infile, /get_lun
readf, lun, s, format ='(a1000)'    ; longest line seen is 600 long
close, lun
free_lun, lun

if keyword_set(ion) then begin
; search for ionisation potential
ionpot = 0.
ionpot_err = 0.
ind = strpos(s, 'IE (evaluated)')
indc = where(ind ge 0, nind)
if nind gt 0 then begin
   ; there is often a definition of this statement in the footer, first is data
   s1 = s[indc[0]]                                 ; the string to extract data from
   p1 = ind[indc[0]]                               ; the position of the identified string - data to the right to that
   p2 = strpos(s1, '.', p1)                        ; search forwards for a decimal point
   p3 = strpos(s1, '>', p2, /reverse_search)       ; search backwards for start of entry
   ps = p3 + 1                                     ; start position
   p4 = strpos(s1, ' ', p2)                        ; search forwards for end of entry
   pe = p4 
   ionpot = double(strmid(s1, ps, pe-ps))
 
   ; see if there is an error estimate
   p5 = strpos(s1, '&plusmn;', pe)
   if p5 gt 0 then begin
      ps = strpos(s1, ' ', p5) + 1
      pe = strpos(s1, '<', ps+1)
     ionpot_err = double(strmid(s1, ps, pe-ps))
   endif
endif else begin
   ind = strpos(s, 'IE (eV)')
   indc = where(ind ge 0, nind) + 1
   if nind gt 0 then begin
     s1 = s[indc[0]]
     ipstr = extract_col(s1, 1)
     ; see if there is an error estimate
     p5 = strpos(ipstr, '&plusmn;')
     if p5 gt 0 then begin
       ps = p5 + 8
       ionpot_err = double(strmid(ipstr, ps, 10))
       ionpot = double(strmid(ipstr, 0, ps))
     endif else begin
       ionpot = double(s1)
     endelse
   endif else begin
     print, 'WARNING: Ionisation Energy not found'
   endelse
endelse

;print, ionpot, ionpot_err, format = '(2f10.7)'
endif


; now search for molecular states and constants
ind = strpos(s, '<td nowrap="nowrap"')
cnt = strarrcnt(s, '</td>')                 ; counts number of columns
indc = where(ind ge 0 and cnt gt 12, nind)
if nind le 0 then begin
   print, 'No states found in ', infile   
   goto, skip
endif

if keyword_set(file) then openw, 1, file
j = 0
for i = 0, nind-1 do begin
    ; extract label in first column
    s1 = s[indc[i]]
    label[j] = relabel(extract_col(s1, 1))
    Te_str[j] = clean_number(extract_col(s1, 2))
    w_str[j] = clean_number(extract_col(s1, 3))
    wx_str[j] = clean_number(extract_col(s1, 4))
    wy_str[j] = clean_number(extract_col(s1, 5))
    B_str[j] = clean_number(extract_col(s1, 6))
    alf_str[j] = clean_number(extract_col(s1, 7))
    gam_str[j] = clean_number(extract_col(s1, 8))
    D_str[j] = clean_number(extract_col(s1, 9))
    bet_str[j] = clean_number(extract_col(s1, 10))
    re_str[j] = clean_number(extract_col(s1, 11))
    
    Te_ref[j] = 'HH'
    we_ref[j] = 'HH'
    wxe_ref[j] = 'HH'
    wye_ref[j] = 'HH'
    Be_ref[j] = 'HH'
    alfe_ref[j] = 'HH'
    game_ref[j] = 'HH'
    De_ref[j] = 'HH'
    bete_ref[j] = 'HH'
    re_ref[j] = 'HH'
    
    get_symmetry, label[j], name, lambda, spin, plusmin, oddeven
    nm[j] = name
    lmbda[j] = lambda
    spn[j] = spin
    plsmn[j] = plusmin
    evenodd[j] = oddeven
    

; check if there are any states with same label in lines below
    
    if i eq nind-1 then begin
       jmax = indc[i] + 10   ; search next ten lines if last label
    endif else begin
       jmax = indc[i+1]-1
    endelse
    ss = s[indc[i]+1: jmax]
    cnt2 = strarrcnt(ss, '</td>')
    indc2 = where(cnt2 gt 11, nind2)
    
    if nind2 gt 0 then begin                 ; same state found - note this means multiplicities should go to 1
        spn[j] = 1
    endif else begin
        spn[j] = spin
    endelse
    if keyword_set(file) then printf, 1, label[j], name, lambda, spn[j], plusmin, oddeven, Te_str[j], w_str[j], wx_str[j], wy_str[j], B_str[j], alf_str[j], gam_str[j], D_str[j], bet_str[j], re_str[j], format = '(a40, a7, i3, i3, a3, a3, 10a12)'
    j = j+1    ; increment j
    
    
    if nind2 gt 0 then begin                 ; same state found
       for jj = 0, nind2-1 do begin
       s2 = ss[indc2[jj]]
       if (clean_number(extract_col(s2, 1)) ne '.' or clean_number(extract_col(s2, 2)) ne '.') then begin
       
         print, 'found an extra state in '+infile+'  : '+clean_number(extract_col(s2, 1))
          
          label[j] = relabel(extract_col(s1, 1))
          Te_str[j] = clean_number(extract_col(s2, 1))
          if clean_number(extract_col(s2, 1)) eq '.' then Te_str[j] = Te_str[j-1]
          w_str[j] = clean_number(extract_col(s2, 2))
          wx_str[j] = clean_number(extract_col(s2, 3))
          wy_str[j] = clean_number(extract_col(s2, 4))
          B_str[j] = clean_number(extract_col(s2, 5))
          alf_str[j] = clean_number(extract_col(s2, 6))
          gam_str[j] = clean_number(extract_col(s2, 7))
          D_str[j] = clean_number(extract_col(s2, 8))
          bet_str[j] = clean_number(extract_col(s2, 9))
          re_str[j] = clean_number(extract_col(s2, 10))
    
          Te_ref[j] = 'HH'
          we_ref[j] = 'HH'
          wxe_ref[j] = 'HH'
          wye_ref[j] = 'HH'
          Be_ref[j] = 'HH'
          alfe_ref[j] = 'HH'
          game_ref[j] = 'HH'
          De_ref[j] = 'HH'
          bete_ref[j] = 'HH'
          re_ref[j] = 'HH'
    
          get_symmetry, label[j], name, lambda, spin, plusmin, oddeven
          nm[j] = name
          lmbda[j] = lambda
          spn[j] = 1 ; spin always 1 if term is divided up
          plsmn[j] = plusmin
          evenodd[j] = oddeven
          if keyword_set(file) then printf, 1, label[j], name, lambda, spn[j], plusmin, oddeven, Te_str[j], w_str[j], wx_str[j], wy_str[j], B_str[j], alf_str[j], gam_str[j], D_str[j], bet_str[j], re_str[j], format = '(a40, a7, i3, i3, a3, a3, 10a12)'
    
          j = j+1
        endif
        endfor
    endif


endfor
if keyword_set(file) then close, 1  

; check through for states with same name side by side and set 2S+1 to 1 in these cases
; the state is divided into its components

for i = 0, nstr - 2 do begin
   if  (nm[i+1] eq nm[i]) then begin
    spn[i] = 1
    spn[i+1] = 1
   endif
endfor


skip:
  
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

numstruct = { $
    label  : label,       $
    name   : nm,          $
    lambda : lmbda,       $
    spin   : spn,         $
    plusmin: plsmn,       $
    oddeven: evenodd,       $
    Te   :    double(clean_Te(Te_str)),      $
    we   :    double(clean_num(w_str)),        $
    wxe  :    double(clean_num(wx_str)),       $
    wye  :    double(clean_num(wy_str)),       $
    Be   :    double(clean_num(B_str)),        $
    alfe :    double(clean_num(alf_str)),      $
    game :    double(clean_num(gam_str)),      $
    De   :    double(clean_num(D_str)),        $
    bete :    double(clean_num(bet_str)),      $
    re   :    double(clean_num(re_str))         $
}

end 
