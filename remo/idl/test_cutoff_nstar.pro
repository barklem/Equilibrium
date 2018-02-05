; script: test calculation of partition function


; element/ion symbol arrays
elemSymArr = ['Fe', 'C', 'O', 'Si', 'Na', 'Mg', 'Mg', 'Ca', 'Ca', 'Ba', 'Y',  'Li', 'Al', 'K', 'Rb', 'Cs', 'Sr', 'Ba']
ionSymArr  = ['I',  'I', 'I', 'I',  'I',  'I',  'II', 'I',  'II', 'II', 'II', 'I',  'I',  'I', 'I',  'I',  'I',  'I' ]
nelem = n_elements(elemSymArr)

; cutoff effective principal quantum number
nstar_cutoff_arr = [ 10., 12., 15., 17., 20.]
ncutoff = n_elements(nstar_cutoff_arr)


; directories
main_dir = '~/Dropbox/Shared/Equilibrium'
atom_dir = main_dir+'/remo/molcalc_test'
plots_dir = main_dir+'/remo/plots'

; some physical constants
hplanck = 4.135667662d-15        ; eV s
kboltz  = 8.6173303d-5           ; eV K-1
cc = 2.99792458d10               ; cm s-1
; ionisation energy HI
eion_HI = 13.6                   ; eV
; Rydberg constant
rydberg_inf = hplanck*cc/eion_HI ; cm-1

; temperature array
temp = [1e-5,  1e-4,  1e-3,  1e-2,  0.1,   0.15,  0.2,   0.3,   0.5,   0.7, $
        1.0,   1.3,   1.7,   2.0,   3.0,   5.0,   7.,    10.,   15.,   20., $
        30.,   50.,   70.,   100.,  130.,  170.,  200.,  250.,  300.,  500., $
        700.,  1000., 1500., 2000., 3000., 4000., 5000., 6000., 7000., 8000., $
        9000., 10000. ]
ntemp = n_elements(temp)


; plot PS parameters
encapsulated  = 1
landscape     = 0
font          = -1
default_thick = 4
charsize      = 1.4
charthick     = default_thick
position      = [0.15,0.15,0.9,0.9] ;position main window
position_inset= [0.52,0.25,0.87,0.53] ;position of inset

; plot labels
xtitle = 'Temperature [K]'
ytitle = 'partition function ratio (partial/total)'

; plot range
xmin = 1.0e2
xmax = 1.2e4
ymin = 0.0
ymax = 1.2

; plot, other parameters
thick = 6
colorStrArr = [ 'blue', 'blu3', 'tomato', 'grn4', 'medium orchid', 'gold', 'turquoise', 'green', 'deep pink' ]


; main
cd, current=current_dir

for ielem=0L,nelem-1 do begin
;for ielem=1,1 do begin

; element/ion symbol
   elemSym = elemSymArr[ielem] 
   ionSym = ionSymArr[ielem]
   elemIon = elemSym+'_'+ionSym

; plot title
   title = elemSym+' '+ionSym

; input element/ion file
; output PS file
   atom_file = elemIon+'_nist.atom'
   ps_file = elemIon+'_partf_nstar_cutoff.eps'

; read element/ion data
   cd, atom_dir
   read_atom, atom_file, jdata, wavenum, eion
   ndata = n_elements(wavenum)
   edata = hplanck*cc*wavenum
   gdata = 2*jdata + 1

; effective nuclear charge ( "Z-N" )
   case ionSym of
      'I'  : nuclear_chrg_eff = 1
      'II' : nuclear_chrg_eff = 2
      'III': nuclear_chrg_eff = 3      
   endcase

; select levels with energies below eion
   ww = where(edata lt eion, nww)

; effective principal quantum number
   nstar = nuclear_chrg_eff / sqrt((eion - edata[ww])/eion_HI)

; compute partition function, "total"
   partf = fltarr(ntemp)
   for j=0L,nww-1 do begin
      k=ww[j]
      partf = partf + gdata[k] * exp( -edata[k]/kboltz/temp )
   endfor

; compute partition function, "partial" (cutoff)
   partf_alt = fltarr(ntemp, ncutoff)
   for icutoff=0L,ncutoff-1 do begin
; where effective principal quantum number is less than nstar_cutoff
      nstar_cutoff = nstar_cutoff_arr[icutoff]
      ww_cutoff = where(nstar le nstar_cutoff, nww_cutoff)
      for i=0L,nww_cutoff-1 do begin
         k = ww[ww_cutoff[i]]
         partf_alt[*,icutoff] = partf_alt[*,icutoff] + gdata[k] * exp( -edata[k]/kboltz/temp )
      endfor
   endfor

; open plot
   cd, plots_dir
   cgPS_open, ps_file, encapsulated=encapsulated, landscape=landscape, $
              font=font, charsize=charsize, default_thick=default_thick, $
              /nomatch

; plot frame
   cgplot, temp, fltarr(ntemp)+1.0, $
           title=title, xtitle=xtitle, ytitle=ytitle, $
           XRange=[xmin,xmax], YRange=[ymin,ymax], /Xlog, $
           linest=2, thick=thick, color='gray'

; plot partition function ratios
   for icutoff=0L,ncutoff-1 do begin
      colorStr = colorStrArr[icutoff]
      cgPlot,/ov,temp,partf_alt[*,icutoff]/partf, $
             color=colorStr, thick=thick
   endfor

; legend
   dxlocl = (!x.crange[1]-!x.crange[0])*0.05
   dyloc  = (!y.crange[1]-!y.crange[0])*0.06

   xloc0 = 10.^(!x.crange[0]+dxlocl)
   yloc0 = 0.6
   titles=string(nstar_cutoff_arr,format='(I02)')
   cgtext,xloc0,yloc0+dyloc,'n* cutoff'
   cglegend, titles=titles, color=colorStrarr[0:ncutoff-1], lines=lonarr(ncutoff), $
             location = [xloc0,yloc0], $
             /data, charthick=charthick, charsize=charsize, thick=thick

   if (nstar[nww-1] ge max(nstar_cutoff_arr)) then color_str_bot='black' else color_str_bot='tomato'
   cgtext,xloc0,0.05,'(max n* in data ='+string(nstar[nww-1],format='(I3)')+')',color=color_str_bot

; plot inset, zoom-in
   cgplot,/NoErase,temp,partf_alt[*,0]/partf, $
          position=position_inset, $
          thick=thick, $
          XRange=[4000,7000],YRange=[0.89,1.01], $
          XMinor=2,YMinor=2,XTicks=3,XStyle=2
; plot partition function ratios
   for icutoff=0L,ncutoff-1 do begin
      colorStr = colorStrArr[icutoff]
      cgPlot,/ov,temp,partf_alt[*,icutoff]/partf, $
             position=position_inset, $
             color=colorStr, thick=thick
   endfor

   
; close plot
   cgPS_close

; convert EPS to PDF
spawn,'epstopdf '+ps_file

; end loop over elements
endfor

;end
cd, current_dir
end
