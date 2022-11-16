; Elements and ionization stages, symbols
elem_sym = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', $
               'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', $
               'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', $
               'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', $
               'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', $
               'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', $
               'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', $
               'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', $
               'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', $
               'Pa', 'U' ]
nelem    = n_elements( elem_sym )

ion_sym  = ['I', 'II', 'III']
nion     = n_elements( ion_sym )

for i=0L,nelem-1 do begin
   
   for j=0L,nion-1 do begin
      elem_ion = string( elem_sym[i], '_', ion_sym[j], format='(A2,A1,A-3)')
      fname    = string( elem_sym[i], '_', ion_sym[j], format='(A2,A1,A)') + '_nist.atom'
      print, format='(A7,4X,A)', elem_ion, fname
   endfor
endfor
end
