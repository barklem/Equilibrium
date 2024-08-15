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

end

