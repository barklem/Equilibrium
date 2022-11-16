pro init_atom_elem

@common_atom_elem

; atomic elements: symbols
  elemSymbols = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', $
                 'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', $
                 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', $
                 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', $
                 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', $
                 'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', $
                 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', $
                 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', $
                 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', $
                 'Pa', 'U' ]

; 92 atomic elements ...
  numElem     = size( elemSymbols, /n_elements)

; atomic numbers
  atomNum  = lindgen(numElem)+1L
  
  
; atomic masses (amu) for the first 92 elements; data from NIST:
; http://physics.nist.gov/cgi-bin/Compositions/
;        stand_alone.pl?ele=&all=all&ascii=ascii&isotype=some
;
  atomMass =  [ $
              1.00794,   4.002602,  6.941,     9.012182,   10.811,    $
              12.0107,   14.0067,   15.9994,   18.9984032, 20.1797,   $
              22.98977,  24.3050,   26.981538, 28.0855,    30.973761, $
              32.065,    35.453,    39.948,    39.0983,    40.078,    $
              44.95591,  47.867,    50.9415,   51.9961,    54.938049, $
              55.845,    58.93320,  58.6934,   63.546,     65.409,    $
              69.723,    72.64,     74.92160,  78.96,      79.904,    $
              83.798,    85.4678,   87.62,     88.90585,   91.224,    $
              92.90638,  95.94,     98.,       101.07,     102.90550, $
              106.42,    107.8682,  112.411,   114.818,    118.710,   $
              121.760,   127.60,    126.90447, 131.293,    132.90545, $
              137.327,   138.9055,  140.116,   140.90765,  144.24,    $
              145,       150.36,    151.964,   157.25,     158.92534, $
              162.500,   164.93032, 167.259,   168.93421,  173.04,    $
              174.967,   178.49,    180.9479,  183.84,     186.207,   $
              190.23,    192.217,   195.078,   196.96655,  200.59,    $
              204.3833,  207.2,     208.98038, 209,        210,       $
              222,       223,       226,       227,        232.0381,  $
              231.03588, 238.02891 ]


; indices of alpha elements: O, Ne, Mg, Si, S, Ar, Ca, Ti
  is_alphaElem    = lonarr(numElem)
  is_alphaElem[*] = 0
  ww_alphaElem    = [ 8, 10, 12, 14, 16, 18, 20, 22 ] -1
  is_alphaElem[ww_alphaElem] = 1

; indices of metals
  is_metalElem      = lonarr(numElem)
  is_metalElem[*]   = 1
  is_metalElem[0:1] = 0
  ww_metalElem      = where(is_metalElem eq 1)

end
