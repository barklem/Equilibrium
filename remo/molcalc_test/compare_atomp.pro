restore,'alldata.sav',/verb
read_sta, aiST, cST, maxE
read_atomp_irwin,ind,code,label,eion,nco,co,label2

; change aiST labels to conform to standard used for label2 and atomid
;
numst = size(aist,/n_elem)
aist2 = strarr(numst)

for i=0,numst-1 do begin
	len = 0
	pos = 0
	el  = ''
	s1  = aist[i]
	pos = stregex(s1,'\++',len=len)
	if (pos ne -1) then begin
	  el  = strmid(s1,0,pos)
	  el  = strupcase(el)+'_I'
	  for j=0,len-1 do begin
	     el = el+'I'
	  endfor
	endif else begin
	  el = strupcase(s1)+'_I'
	endelse	
	aist2[i] = el
endfor


;
; read atomdata
;
;
atomdfile = '/Users/remo/convec/atomdata_orig'

text 	= ''
num_ad	= 92L


label_ad = strarr(num_ad*3L)
eion_ad  = fltarr(num_ad*3L)
qtab_ad	 = fltarr(3,num_ad*3L)
qmin_ad  = fltarr(num_ad*3L)
tmin_ad	 = fltarr(num_ad*3L)

j = 0
openr, lun, atomdfile, /get_lun

   readf, lun, text		;read first line containing comment
   for i=0L,num_ad-1L do begin

      readf, lun, format='(a2,x,f5.2,x,f5.2,x,f5.2,2x,f4.2)',el,chi1,chi2,chi3,iflag	

      el = strcompress(el,/remove_all)

      readf, lun, format='(31x,f5,x,f5,x,f5,x,f5,x,f5)',part1,part2,part3,tmin,part4	
      label_ad[j] = el+'_I'
      eion_ad[j]  = chi1
      qtab_ad[*,j]= [part1,part2,part3]
      tmin_Ad[j]  = tmin
      qmin_ad[j]  = part4

      j++

      readf, lun, format='(31x,f5,x,f5,x,f5,x,f5,x,f5)',part1,part2,part3,tmin,part4

      label_ad[j] = el+'_II'
      eion_ad[j]  = chi2
      qtab_ad[*,j]= [part1,part2,part3]
      tmin_Ad[j]  = tmin
      qmin_ad[j]  = part4

      j++

      readf, lun, format='(31x,f5,x,f5,x,f5,x,f5,x,f5)',part1,part2,part3,tmin,part4

      label_ad[j] = el+'_III'
      eion_ad[j]  = chi3
      qtab_ad[*,j]= [part1,part2,part3]
      tmin_Ad[j]  = tmin
      qmin_ad[j]  = part4

      j++

   endfor
   
free_lun, lun





;
; compare partition function data
;
;

natomid = size(atomid,/n_elem)
theta	= 5040./t

set_plot,'ps'
device,file='~/compare_atom_partf.ps',/color, bit=8
loadct,39

!x.thick=3
!y.thick=3
!p.charsize=1.4
!p.charthick=3

for n=0L,natomid-1L do begin
;for n=98L,100L do begin

  ymax = max(qatom[n,*]) * 1.4
  ymin = min(qatom[n,*]) * 0.8
  
  plot,t/1000.,qatom[n,*], 	 	$
  	title=atomid[n],	$
	ytitle='Partition function',	$
	xtitle='T [1000 K]', $
	xr=[1, 16], /xst, yr=[ymin,ymax], $
	co=0, th=6,linest=0
	xyouts, 2.,ymin+(ymax-ymin)*0.9, 'NIST', charthick=3,charsize=1.4
	
  w1 = where(label2 eq strupcase(atomid[n]))
  w2 = where(aist2  eq strupcase(atomid[n]))
  w3 = where(label_ad  eq strupcase(atomid[n]))
  print, atomid[n],n,w1[0],w2[0],w3[0]
  if (w1 ne -1) then begin
  	oplot, t/1000., exp(poly(alog(t),co[*,w1])), th=6, co=100, linest=2
	xyouts, 2.,ymin+(ymax-ymin)*0.8, 'IRWIN 1981', charthick=3,charsize=1.4,col=100
  endif
  if (w2 ne -1) then begin
  	oplot, t/1000., 10.^(poly(alog10(theta),cst[*,w2])),th=6, co=200, linest=3
	xyouts, 2.,ymin+(ymax-ymin)*0.7, 'SAUVAL & TATUM 1984', charthick=3,charsize=1.4,col=200
  endif
  if (w3 ne -1) then begin
  	temp = [3600., 5700., 8000.]
	qtab = alog( reform( qtab_ad[*,w3] ) )
	qmin = reform( qmin_ad[w3] )
	tmin = tmin_ad[w3] 
	f11 = (qtab[1]-qtab[0])/( temp[1]-temp[0] ) 
	f12 = (qtab[2]-qtab[1])/( temp[2]-temp[1] )
	f2  = (f12-f11)/(temp[2]-temp[0])
	qout = qtab[0] + f11* (t-temp[0]) + f2* (t-temp[0])*(t-temp[1])
	qout = exp(qout)
	ww   = where(t lt tmin[0])
	qout[ww] = qmin
	help,ww
  	oplot, t/1000., qout,th=6, co=150, linest=4
	oplot, temp/1000, exp(qtab), ps=8,syms=1.5,co=150
	xyouts, 2.,ymin+(ymax-ymin)*0.6, 'ATOMDATA (ORIGINAL)', charthick=3,charsize=1.4,col=150
  endif
  
 ; wait,0.5
endfor

!p.thick=1
!x.thick=1
!y.thick=1
!p.charsize=1
!p.charthick=1
device, /close
set_plot,'x'

end
