pro read_diss, name, adopted, hh=hh, luo=luo, g2=g2, other=other, components=components

file = 'Dissociation/diss.txt'

openr, 1, file

n = file_line(file) 

hh = dblarr(n,2)
luo = dblarr(n,2)
g2 = dblarr(n,2)
other = dblarr(n,2)
adopted = dblarr(n,2)
name = strarr(n)
components = strarr(n,2)

i=0
for j= 0, n-1 do begin
   s = ' '
   readf, 1, s, format = '(a200)'
   s1 = strsplit(s, ' ', /extract)
   name(i) = s1(0)
   if s1(0) eq '#' then begin
     n = n-1
   endif else begin   
     components(i,0)=s1(1)
     components(i,1)=s1(2)
     if strpos(s1[3], '-') ge 0 then hh(i,0) = 0.d0 else hh(i,0) = double(s1[3])
     if strpos(s1[4], '-') ge 0 then hh(i,1) = 0.d0 else hh(i,1) = double(s1[4])
     if strpos(s1[5], '-') ge 0 then luo(i,0) = 0.d0 else luo(i,0) = double(s1[5])
     if strpos(s1[6], '-') ge 0 then luo(i,1) = 0.d0 else luo(i,1) = double(s1[6])
     if strpos(s1[8], '-') ge 0 then g2(i,0) = 0.d0 else g2(i,0) = double(s1[8])
     if strpos(s1[9], '-') ge 0 then other(i,0) = 0.d0 else other(i,0) = double(s1[9])
     case s1[10] of
        'H': adopted(i,*) = hh(i,*)
        'L': adopted(i,*) = luo(i,*)
        'C': adopted(i,*) = g2(i,*)
        'O': adopted(i,*) = other(i,*)
     endcase
     i = i+1
   endelse
endfor

name = name(0:n-1)
components = components(0:n-1,*)
hh = hh(0:n-1,*)
luo = luo(0:n-1,*)
g2 = g2(0:n-1,*)
other = other(0:n-1,*)
adopted = adopted(0:n-1,*)

close, 1
end