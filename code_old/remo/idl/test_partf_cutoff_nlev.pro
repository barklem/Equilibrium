; script: test calculation of partition function


; element/ion symbol arrays
elemSymArr = ['Fe', 'C', 'O', 'Si', 'Na', 'Ca', 'Ca', 'Ba', 'Y']
ionSymArr  = ['I',  'I', 'I', 'I',  'I',  'I',  'II', 'II', 'II' ]
nelem = n_elements(elemSymArr)


; cutoff levels: fraction of total number of levels
cutoffArr = [0.95, 0.9, 0.75, 0.5, 0.25, 0.1]
ncutoff = n_elements(cutoffArr) 

; directories
main_dir = '~/Dropbox/Shared/Equilibrium'
atom_dir = main_dir+'/remo/molcalc_test'
plots_dir = main_dir+'/remo/plots'

; some physical constants
hplanck = 4.135667662d-15       ; eV s
kboltz  = 8.6173303d-5          ; eV K-1
cc = 2.99792458d10              ; cm s-1

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
colorStrArr = [ 'blu3', 'tomato', 'grn4', 'medium orchid', 'gold', 'turquoise', 'green', 'deep pink' ]


;;elemIonArr = ['Fe_I']
;;nelemIonArr = n_elements(elemIonArr)


; main
cd, current=current_dir

for ielem=0L,nelem-1 do begin

; element/ion symbol
   elemSym = elemSymArr[ielem] 
   ionSym = ionSymArr[ielem]
   elemIon = elemSym+'_'+ionSym

; plot title
   title = elemSym+' '+ionSym

; input element/ion file
; output PS file
   atom_file = elemIon+'_nist.atom'
   ps_file = elemIon+'_partf_cutoff_nlev.eps'

; read element/ion data
   cd, atom_dir
   read_atom, atom_file, jdata, wavenum, eion
   ndata = n_elements(wavenum)
   edata = hplanck*cc*wavenum
   gdata = 2*jdata + 1

; compute partition function
   partf = fltarr(ntemp)
   ww = where(edata le eion, nww)
   nlevtot = nww
   for j=0L,nlevtot-1 do begin
      k=ww[j]
      partf = partf + gdata[k] * exp( -edata[k]/kboltz/temp )
   endfor
   

; compute modified partition functions, assuming a cutoff excitation
; energy
   partf_alt = fltarr(ntemp,ncutoff)

   for icutoff=0L,ncutoff-1 do begin
      cutoff = cutoffArr[icutoff]
      ecutoff = cutoff * eion
      nlev = long(nlevtot * cutoff)
      for j=0L,nlev-1 do begin
         k=ww[j]
         partf_alt[*,icutoff] = partf_alt[*,icutoff] + gdata[k] * exp( -edata[k]/kboltz/temp )
      endfor
   endfor

; plot
; open PS
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
   yloc0 = 0.5
   xloc1 = 10.^(!x.crange[0]+dxlocl*4.)
   yloc1 = yloc0 - dyloc
   cgText,xloc0,yloc0,'Cutoff (fraction of total number of levels)'
   for icutoff=0L,ncutoff-1 do begin
      xloc = xloc1
      yloc = yloc1 - icutoff * dyloc
      cutoff = cutoffArr[icutoff]
      cgText, xloc, yloc, string(cutoff*100.,format='(I3)')+' %'
      colorStr = colorStrArr[icutoff]
      cgPlot,/ov,10.^(!x.crange[0]+dxlocl*[1.,3.5]),[1,1]*(yloc+dyloc*0.3),linest=0,thick=thick,color=colorStr
   endfor

; close
   cgPS_close

endfor

;end
cd, current_dir
end
