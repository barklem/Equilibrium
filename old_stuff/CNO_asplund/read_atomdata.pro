FUNCTION READ_ATOMDATA, atomfile=atomfile

default,atomfile,'~/conv/tables/atomdata_nist'

openr, lun, atomfile, /get_lun

a=''
tmp=fltarr(3)
tmp2=dblarr(6)

readf, lun, a
readf, lun, natom, nmol

label=strarr(natom)
ion=fltarr(natom,2)
weight=fltarr(natom)
partcoeff=dblarr(natom,3,6)

for i=0, natom-1 do begin
    readf, lun, a, tmp, format='(A2,3F8.3)'
    label[i]=strcompress(a,/remove)
    ion[i,*]=tmp[0:1]
    weight[i]=tmp[2]
    for j=0,2 do begin 
        readf, lun, tmp2
        partcoeff[i,j,*]=tmp2
    endfor
endfor

tmp=fltarr(3)
tmp2=dblarr(6)
tmp3=indgen(4)

mollabel=strarr(nmol)
moldiss=fltarr(nmol)
molpartcoeff=dblarr(nmol,6)
moldefit=dblarr(nmol,6)
moleefit=dblarr(nmol,6)

for i=0, nmol-1 do begin
    readf, lun, a, tmp, tmp3, format='(A8,3F8.3,4I3)'
    mollabel[i]=strcompress(a,/remove)
    moldiss[i,*]=tmp[0]
    readf, lun, tmp2
    molpartcoeff[i,*]=tmp2
    if tmp3[0] gt 0.0 then begin 
       tmp4=fltarr(tmp3[2])
       readf, lun, tmp4
       moldefit[i,0:tmp3[2]-1]=tmp4
    endif
    if tmp3[1] gt 0.0 then begin
       tmp5=fltarr(tmp3[3])
       readf, lun, tmp5
       moleefit[i,0:tmp3[3]-1]=tmp5
    endif
endfor

free_lun, lun

data={natom:natom,nmol:nmol,label:label,ion:ion,weight:weight,partcoeff:partcoeff,$
      mollabel:mollabel,moldiss:moldiss,molpartcoeff:molpartcoeff,moldefit:moldefit,moleefit:moleefit}

return, data

END


PRO PLOT_ATOMDATA, atomfile=atomfile, ps=ps

default,atomfile,'~/conv/tables/atomdata_nist'

if keyword_set(ps) then pscol, 'Q.ps', xsize=18, ysize=24, xoff=1, yoff=1

@basic_colors
!p.multi=[0,1,3]

data=read_atomdata(atomfile=atomfile)

temp=findgen(15)*1000.+1000.
lnt=alog(temp)
theta=5039.75/temp

logq=fltarr(data.natom,3,n_elements(temp))
q=logq
sion=['I','II','III']

for i=0,data.natom-1 do begin
;for i=4,8 do begin
    for j=0,2 do begin
        ;logq[i,j,*]=data.partcoeff[i,j,0]+data.partcoeff[i,j,1]*theta^1.+ $
        ;          data.partcoeff[i,j,2]*theta^2.+data.partcoeff[i,j,3]*theta^3.+$
        ;          data.partcoeff[i,j,4]*theta^4.+data.partcoeff[i,j,5]*theta^5.
        ;q[i,j,*]=10.^logq[i,j,*]
        logq[i,j,*]=data.partcoeff[i,j,0]+data.partcoeff[i,j,1]*lnt^1.+ $
                  data.partcoeff[i,j,2]*lnt^2.+data.partcoeff[i,j,3]*lnt^3.+$
                  data.partcoeff[i,j,4]*lnt^4.+data.partcoeff[i,j,5]*lnt^5.
        q[i,j,*]=exp(logq[i,j,*])
        ;help, q, theta
        plot, temp, q[i,j,*], /nodata, title=data.label[i]+sion[j], charsize=1.6
        oplot, temp, q[i,j,*], col=red
        ind=where(temp eq 3000.)
        xyouts, temp[1], !y.crange[0]*0.1+!y.crange[1]*0.9, 'U(T=3000) ='+string(q[i,j,ind], format='(F8.3)')
        ind=where(temp eq 5000.)
        xyouts, temp[1], !y.crange[0]*0.3+!y.crange[1]*0.7, 'U(T=5000) ='+string(q[i,j,ind], format='(F8.3)')
        ind=where(temp eq 8000.)
        xyouts, temp[1], !y.crange[0]*0.5+!y.crange[1]*0.5, 'U(T=8000) ='+string(q[i,j,ind], format='(F8.3)')
        ind=where(temp eq 12000.)
        xyouts, temp[1], !y.crange[0]*0.7+!y.crange[1]*0.3, 'U(T=12000)='+string(q[i,j,ind], format='(F8.3)')
    endfor
    ;wait, 1
endfor

if keyword_set(ps) then psend

END

FUNCTION PARTITIONFUNCTION, temp, z, ion, atomfile=atomfile, eion=eion, label=label

default,atomfile,'~/conv/tables/atomdata_nist'

data=read_atomdata(atomfile=atomfile)

lnt=alog(temp)
logq=data.partcoeff[z-1,ion,0]+data.partcoeff[z-1,ion,1]*lnt^1.+ $
     data.partcoeff[z-1,ion,2]*lnt^2.+data.partcoeff[z-1,ion,3]*lnt^3.+$
     data.partcoeff[z-1,ion,4]*lnt^4.+data.partcoeff[z-1,ion,5]*lnt^5.
q=exp(logq)

;print, data.label[z-1], reform(data.ion[z-1,*]) 
eion=data.ion[z-1,ion]
label=data.label[z-1]

return, q

END

PRO WRITE_PARTITIONFUNCTION

temp=[3000.,5000.,8000.,12000.]
ion=['I', 'II', 'III']
openw, gun, '~/latex/articles/sun_araa/partition_functions/partition.tex', /get_lun

for i=3, 92 do begin
   for j=0, 1 do begin
      
      part=partitionfunction(temp, i, j, eion=eion, label=label)
      printf, gun, label, ' &', ion[j], ' &',  eion, ' &', part[0], ' &', part[1], ' &', part[2], ' &', part[3], ' \\', $
             format='(A4,A2,A4,A2,F8.3,A2,F10.2,A2,F10.2,A2,F10.2,A2,F10.2,A3)'
      print, label, ' &', ion[j], ' &',  eion, ' &', part[0], ' &', part[1], ' &', part[2], ' &', part[3], ' \\', $
             format='(A4,A2,A4,A2,F8.3,A2,F10.2,A2,F10.2,A2,F10.2,A2,F10.2,A3)'
   endfor
endfor

free_lun, gun

END