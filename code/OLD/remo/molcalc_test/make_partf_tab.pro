;
; read ionization potentials for the first 92 elements, three ionization stages
;

dir      = '~/People/barklem/Equilibrium/remo/molcalc_test'
eionfile = 'ionpot_moog_all.txt'
alldatafile = 'alldata.sav'

readcol, dir+'/'+eionfile, atomnum, el, eion1, eion2, eion3, format='i,a,f,f,f', comment='*'
nel = size(atomnum,/n_elem)


;
; define and initialize table with coefficients of the polynomial fits (all elements, three ionization stages)
;

ncoef	  = 6L
tabco	  = dblarr(ncoef, nel*3L)
tabco[*,*]= 0.0d0
label_ion = strarr(nel*3L)

ntemp	= 551
temp	= dindgen(ntemp)/double(ntemp-1)*15000.+1000.
test1   = dindgen(ntemp,3L*nel)
test2   = test1
test3   = test1


;
; atomic mass (amu) for the first 92 elements; from NIST http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii&isotype=some
;

atom_mass =  [ $
	1.00794,   4.002602,  6.941,     9.012182, 10.811,    12.0107, 14.0067,   15.9994, 18.9984032, 20.1797, $
	22.989770, 24.3050,   26.981538, 28.0855,  30.973761, 32.065,  35.453,    39.948,  39.0983,    40.078,  $
	44.955910, 47.867,    50.9415,   51.9961,  54.938049, 55.845,  58.933200, 58.6934, 63.546,     65.409, $
	69.723,    72.64,     74.92160,  78.96,    79.904,    83.798,  85.4678,   87.62,   88.90585,   91.224, $
	92.90638,  95.94,     98,         101.07,   102.90550, 106.42,  107.8682,  112.411, 114.818,    118.710, $
	121.760,   127.60,    126.90447, 131.293,  132.90545, 137.327, 138.9055,  140.116, 140.90765,  144.24, $
	145,       150.36,    151.964,   157.25,   158.92534, 162.500, 164.93032, 167.259, 168.93421,  173.04, $
	174.967,   178.49,    180.9479,  183.84,   186.207,   190.23,  192.217,   195.078, 196.96655,  200.59, $
	204.3833,  207.2,     208.98038, 209,      210,       222,     223,       226,     227,        232.0381, $
	231.03588, 238.02891 ]


; 
; Up-to-date MOOG 2007 partition functions for some of the rare earth elements
;

nmoog		= 27L
label_moog	= strarr(nmoog)
co_moog		= dblarr(6,nmoog)
i=0  &  label_moog[i]='CR_II'  & co_moog[*,i]=   [	2.66384894d+3, -1.68227473d+3,   4.22830259d+2, -5.28071678d+1,  3.27472214d+0, -8.05963460d-2]
i=1  &  label_moog[i]='RB_I'   & co_moog[*,i]=   [	5.03676283d+3, -3.16069504d+3,   7.89175381d+2, -9.79501783d+1,  6.04047992d+0, -1.47981323d-1]
i=2  &  label_moog[i]='MO_II'  & co_moog[*,i]=   [	4.16034835d+3, -2.65480884d+3,   6.74631372d+2, -8.52665787d+1,  5.35745654d+0, -1.33793443d-1]
i=3  &  label_moog[i]='PR_I'   & co_moog[*,i]=   [  -3.04466224d+3,  1.88966140d+3,  -4.64804826d+2,  5.66502603d+1, -3.42059183d+0,  8.19184298d-2   ]
i=4  &  label_moog[i]='PR_II'  & co_moog[*,i]=   [  -1.73962063d+3,  1.07851006d+3,  -2.64849957d+2,  3.22161314d+1, -1.93975232d+0,  4.62805386d-2   ]
i=5  &  label_moog[i]='ND_I'   & co_moog[*,i]=   [  -1.64889203d+3,  9.47177848d+2,  -2.12837179d+2,  2.33051095d+1, -1.23575686d+0,  2.52043511d-2   ]
i=6  &  label_moog[i]='ND_II'  & co_moog[*,i]=   [  -2.06698421d+2,  9.79972763d+1,  -1.55582617d+1,  7.28892861d-1,  3.56166693d-2, -2.95215536d-3   ]
i=7  &  label_moog[i]='GD_I'   & co_moog[*,i]=   [	1.83011048d+3, -1.19539486d+3,   3.10632505d+2, -4.00859062d+1,  2.56842403d+0, -6.53128117d-2]
i=8  &  label_moog[i]='TB_I'   & co_moog[*,i]=   [	5.85584746d+2, -3.89743229d+2,   1.03759088d+2, -1.37350828d+1,  9.03674535d-1, -2.35845042d-2]
i=9  &  label_moog[i]='TB_II'  & co_moog[*,i]=   [  -9.48002791d+2,  6.16808235d+2,  -1.58524150d+2,  2.01610509d+1, -1.26856460d+0,  3.16355813d-2   ]
i=10 &  label_moog[i]='DY_I'   & co_moog[*,i]=   [    3.53878819d+3, -2.31997287d+3,   6.06023177d+2, -7.87382964d+1,  5.08456496d+0, -1.30426644d-1  ]
i=11 &  label_moog[i]='DY_II'  & co_moog[*,i]=   [    2.77061687d+3, -1.78523018d+3,   4.58361896d+2, -5.85341812d+1,  3.71600397d+0, -9.37402524d-2  ]
i=12 &  label_moog[i]='HO_II'  & co_moog[*,i]=   [    1.34187596d+3, -8.52533583d+2,   2.15784324d+2, -2.71271661d+1,  1.69281185d+0, -4.19020326d-2  ]
i=13 &  label_moog[i]='ER_I'   & co_moog[*,i]=   [    2.02411710d+3, -1.35836078d+3,   3.62878018d+2, -4.81385931d+1,  3.16795944d+0, -8.26308027d-2  ]
i=14 &  label_moog[i]='ER_II'  & co_moog[*,i]=   [    1.62058945d+3, -1.06472504d+3,   2.78625547d+2, -3.62193240d+1,  2.33715989d+0, -5.98219225d-2  ]
i=15 &  label_moog[i]='YB_I'   & co_moog[*,i]=   [    7.47839800d+3, -4.67427018d+3,   1.16206341d+3, -1.43581759d+2,  8.81293038d+0, -2.14850223d-1  ]
i=16 &  label_moog[i]='LU_II'  & co_moog[*,i]=   [    2.44026764d+3, -1.59016793d+3,   4.12298445d+2, -5.31395892d+1,  3.40232351d+0, -8.65004469d-2  ]
i=17 &  label_moog[i]='TH_I'   & co_moog[*,i]=   [    4.76054139d+2, -3.25783548d+2,   8.97134083d+1, -1.23483587d+1,  8.46088455d-1, -2.29570332d-2  ]
i=18 &  label_moog[i]='TH_II'  & co_moog[*,i]=   [  -1.22179869d+3,  7.59698688d+2,  -1.86944859d+2,  2.27591989d+1, -1.37021999d+0,  3.26865578d-2   ]
i=19 &  label_moog[i]='U_II'   & co_moog[*,i]=   [  -2.76709404d+2,  1.30426239d+2,  -2.06221054d+1,  9.83307627d-1,  4.19297328d-2, -3.60428836d-3   ]
i=20 &  label_moog[i]='EU_III' & co_moog[*,i]=   [  -2.02846263d+3,  1.24521575d+3,  -3.04725159d+2,  3.72010364d+1, -2.26584074d+0,  5.50893177d-2   ]
i=21 &  label_moog[i]='TB_III' & co_moog[*,i]=   [    8.11558123d+2, -4.67434582d+2,   1.07583056d+2, -1.23251149d+1,  7.02860514d-1, -1.59394115d-2  ]
i=22 &  label_moog[i]='DY_III' & co_moog[*,i]=   [    1.30822682d+3, -7.46679156d+2,   1.69716518d+2, -1.91592590d+1,  1.07376244d+0, -2.38653520d-2  ]
i=23 &  label_moog[i]='HO_III' & co_moog[*,i]=   [    6.98122851d+3, -4.11418184d+3,   9.66307217d+2, -1.13063284d+2,  6.59141456d+0, -1.53163679d-1  ]
i=24 &  label_moog[i]='ER_III' & co_moog[*,i]=   [    2.65562213d+3, -1.49666461d+3,   3.34740125d+2, -3.70678109d+1,  2.02963822d+0, -4.38640816d-2  ]
i=25 &  label_moog[i]='TH_III' & co_moog[*,i]=   [    7.10076650d+3, -4.15762991d+3,   9.69365844d+2, -1.12510812d+2,  6.50161887d+0, -1.49626283d-1  ]
i=26 &  label_moog[i]='U_III'  & co_moog[*,i]=   [    5.18480529d+3, -3.02564444d+3,   7.03396006d+2, -8.14028851d+1,  4.68959696d+0, -1.07563314d-1  ]


;
; read coefficients of polynomial fit to the partition functions (Irwin 1981);
; fill in tabco with Irwin 1981 data
;

read_atomp_irwin,ind,code,label,eion,nco,co,label2

w1=-1L									;reset, just in case...
w2=-1L
w3=-1L

for i=0L,nel-1L do begin
	
	label_ion[3L*i+0]=el[i]+'_I'
	label_ion[3L*i+1]=el[i]+'_II'
	label_ion[3L*i+2]=el[i]+'_III'
	
	w1=where(label2 eq strupcase(el[i])+'_I'  )
	w2=where(label2 eq strupcase(el[i])+'_II' )
	w3=where(label2 eq strupcase(el[i])+'_III')
	tabco[*,3L*i]   = double( co[*,w1[0]] )
	tabco[*,3L*i+1] = double( co[*,w2[0]] )
	test1[*,3L*i]   = exp(poly(alog(temp),tabco[*,3L*i]))
	test1[*,3L*i+1] = exp(poly(alog(temp),tabco[*,3L*i+1]))
	if (w3[0] ne -1) then begin
		tabco[*,3L*i+2] = double( co[*,w3[0]] )		; H_III doesn't exist...
		test1[*,3L*i+2] = exp(poly(alog(temp),tabco[*,3L*i+2]))
	endif
endfor



;
; load partition functions computed with NIST data about energy levels;
; compute polynomial fits;
; update tabco when NIST data are available
;

restore, dir+'/'+alldatafile, /verbose


w1=-1L									;reset, just in case...
w2=-1L
w3=-1L

tmin	= 1000.
tmax	= 16000.
ww	= where(t ge tmin and t le tmax)
logt	= alog(t)

for i=0L,nel-1L do begin
	w1=where(atomid eq el[i]+'_I'  )
	w2=where(atomid eq el[i]+'_II' )
	w3=where(atomid eq el[i]+'_III')
	if (w1[0] ne -1) then begin
		cofit1	= poly_fit(logt[ww],alog(qatom[w1[0],ww]), ncoef-1L,/double)
		tabco[*,3L*i]   = double( reform(cofit1,ncoef,1) )
		test2[*,3L*i]   = exp(poly(alog(temp),tabco[*,3L*i]))
	endif
	if (w2[0] ne -1) then begin
		cofit2	= poly_fit(logt[ww],alog(qatom[w2[0],ww]), ncoef-1L, /double)
		tabco[*,3L*i+1] = double( reform(cofit2,ncoef,1) )
		test2[*,3L*i+1] = exp(poly(alog(temp),tabco[*,3L*i+1]))
	endif
	if (w3[0] ne -1) then begin
		cofit3	= poly_fit(logt[ww],alog(qatom[w3[0],ww]), ncoef-1L,/double)	
		tabco[*,3L*i+2] = double( reform(cofit3,ncoef,1) )
		test2[*,3L*i+2] = exp(poly(alog(temp),tabco[*,3L*i+2]))
	endif
endfor

w0	= -1L

test3[*,*]=0.0d0
for i=0L,nmoog-1L do begin
	w0	= where( strupcase(label_ion) eq strupcase(label_moog[i]) )
	if (w0[0] ne -1) then begin 
		print,w0[0],label_ion[w0[0]]
		tabco[*,w0[0]] = double( co_moog[*,i])
		test3[*,w0[0]] = exp(poly(alog(temp),co_moog[*,i]))
	endif
endfor 

;set_plot,'ps'
;device,file='~/test.ps',/color, bit=8
;loadct,39
;!x.thick=3
;!y.thick=3
;!p.charsize=1.4
;!p.charthick=3
;for i=0,275 do begin
;	ymax=max(test1[*,i])*1.2
;	plot,temp,test1[*,i],yr=[0,ymax],title=label_ion[i]
;	oplot,temp,test2[*,i],co=100
;	oplot,temp,test3[*,i],co=200
;endfor
;!p.thick=1
;!x.thick=1
;!y.thick=1
;!p.charsize=1
;!p.charthick=1
;device, /close
;set_plot,'x'


;
; MOLECULES
;

listmol	= [ 'C2',   'CN',    'CH',    'CO',    'O2',    'OH',    'H2',    'NH',    'N2',    'NO',    'MgH',    'HF' ]
massmol = [24.0214, 26.0174, 13.0186, 28.0101, 31.9988, 17.0073, 2.01588, 15.0146, 28.0134, 30.0061, 23.9977,  20.0063] 
nlistmol = size(listmol, /n_elem)
tabcomol = dblarr(ncoef, nlistmol)


wmol	= -1L
for i=0L,nlistmol-1L do begin
	wmol = where(molid eq listmol[i])
	if (wmol[0] ne -1) then begin
		cofit	= poly_fit(logt[ww],alog(qmol[wmol[0],ww]), ncoef-1L,/double)
		tabcomol[*,i] = double( reform(cofit,ncoef,1) )
	endif
endfor


;
; PRINT TABLE WITH FIT COEFFICIENTS 
;

openw, /get, lun, '~/testatom'
printf,lun, '* NAT   NMOL'
printf,lun, nel, nlistmol
;printf,lun, format='(A4,I2,A6,I2,A7,I2)','NAT=',nel,',NMOL=',nlistmol,',NCOEF=',ncoef
for i=0L, nel-1L do begin
	printf, lun, format='(A-2,2X,F6.3,X,F6.3,2X,F7.3)', strupcase( el[i] ), eion1[i],eion2[i],atom_mass[i]
	printf, lun, format='(31X,6(X,G14.7))', tabco[*,3L*i]
	printf, lun, format='(31X,6(X,G14.7))', tabco[*,3L*i+1]
	printf, lun, format='(31X,6(X,G14.7))', tabco[*,3L*i+2]
;	printf, lun, format='(32X,I2,X,6(X,G14.7))',ncoef, tabco[*,3L*i]
;	printf, lun, format='(32X,I2,X,6(X,G14.7))',ncoef, tabco[*,3L*i+1]
;	printf, lun, format='(32X,I2,X,6(X,G14.7))',ncoef, tabco[*,3L*i+2]
endfor
for i=0L, nlistmol-1L do begin
	read_mol,listmol[i]+'_nist_huber79.mol',ms
	dissoc = ms.d
	printf, lun, format='(A-3,2X,F6.3,X,F6.3,2X,F7.3)', listmol[i], dissoc,0.00,massmol[i]
	printf, lun, format='(31X,6(X,G14.7))', tabcomol[*,i]
;	printf, lun, format='(32X,I2,X,6(X,G14.7))',ncoef, tabcomol[*,i]
endfor

free_lun,lun
END
