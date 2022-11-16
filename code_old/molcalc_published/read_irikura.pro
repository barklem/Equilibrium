pro read_irikura, iristruct

file = 'Irikura/irikura.txt'
n = file_line(file) -1 
openr, 1, file
s = ' '
readf, 1, s, format = '(a200)'   ; one header

name = strarr(n)
we = dblarr(n)
wxe = dblarr(n)
wye = dblarr(n)
Be = dblarr(n)
alfe = dblarr(n)

for i = 0, n-1 do begin
   readf, 1, s, format='(a200)'
   s2 = strsplit(s, /extract)
   name(i) = s2(0)
   we(i) = double(s2(1))
   wxe(i) = double(s2(3))
   if s2(5) ne '.' then wye(i) = double(s2(5))
   Be(i) = double(s2(7))
   alfe(i) = double(s2(9))
endfor
close, 1


iristruct = { $
    n    :    n,         $
    name :    name,      $
    we   :    we,        $
    wxe  :    wxe,       $
    wye  :    wye,       $
    Be   :    Be,        $
    alfe :    alfe      $
  }

end