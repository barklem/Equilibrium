; 1. rewrite ionpot_crc.txt to nicer format; 
; 2. write extra file containing ionization energy for first three
; ions only


dir       = '~/Desktop/ionpot/'
inFile    = dir+'ionpot_CRChandbook_ed91.txt'
outFile1  = dir+'ionpot_CRC_all.txt'
outFile2  = dir+'ionpot_CRC_I-III.txt'


; define some useful constants

numAvog      = 6.02214179d-23
eV_to_J      = 1.60217653d-19
kJpmol_to_eV = eV_to_J / numAvog * 1.0d-3


; define element symbols

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

nElemSym     = size( elemSymbols, /n_elements )



; define roman numerals (for ionization stages)

romNumDec    = [  'I',  'II',  'III',  'IV',  'V',  'VI', 'VII', 'VIII', 'IX', 'X' ]
romNums      = [ romNumDec, 'X'+romNumDec, 'XX'+romNumDec ]      



; open files for input/output

openr, lunIn,   /get_lun, inFile 

openw, lunOut1, /get_lun, outFile1
openw, lunOut2, /get_lun, outFile2



; print headers

newline   = string(10B)         ; define newline character
header    = '# Measured ionization energies (eV) for all elements (and several ionization stages)' + newline + $
            '# Source: CRC Handbook of Chemistry and Physics,  91st edition' + newline + $
            '# ' 

header1   = header + newline + $
            '# '   + $
            string( format='(a3,2x,a13,2x,a3)', 'Z', 'Element', 'Sym') + $
            string( format='(30(2x,a11))', romNums[0:29] )

header2   = header + newline + $
            '# -1 flag indicates missing data' + newline + $
            '# '   + newline + $
            '# '   + $
            string( format='(a3,3x,a3)', 'Z', 'Sym') + $
            string( format='(3(2x,a11))', romNums[0:2] )

printf, lunOut1, header1
printf, lunOut2, header2




; process input file data; 
; reformat it before writing it to output file


dumStr = ''                     ; reset dummy string for input

while ~ eof( lunIn ) do begin

                                ; read dummy string from input file

   readf, lunIn, dumStr

                                ; check whether line is empty

   dumStrCmprs = strcompress( dumStr, /remove_all)
   isStrEmpty  = (dumStrCmprs eq '') 


   if ~isStrEmpty then begin

; reset output strings

      outStr1 = ''
      outStr2 = ''

      splitArr   = strsplit( dumStr, /extract )
      nSplitArr  = size( splitArr, /n_elements )
      
      k          = long( splitArr[0] )
      atomNum    = long( splitArr[1] )
      atomName   = splitArr[2]
      atomSym    = splitArr[3]

      outStr1  = string( format='(2x,i3,2x,a13,2x,a3)', atomNum, atomName, atomSym )
      outStr2  = string( format='(2x,i3,3x,a3)',        atomNum, atomSym )


; number of available measurements of ionization potentials

      nIon       = nSplitArr-4
      isNIonPos  = (nIon  gt 0) 

      nIon2      = min( [3L, nIon] )
      isNIonPos2 = (nIon2 gt 0)


; reset ionization energies (I, II, III) for output file 2

      ionPotArr2 = [ '-1', '-1', '-1' ]


; add ionization potential to output strings

      if isNIonPos then begin
         for i=0L,nIon-1 do begin
            outStr1 = outStr1 + string( format='(2x,a11)', splitArr[4+i] )
         endfor
      endif


; add ionization potential to output string 2

      if isNIonPos2 then begin
         for i=0L,nIon2-1L do begin              
            ionPotArr2[i]  = splitArr[4+i]
         endfor
      endif
      outStr2 = outStr2 + string( format='(3(2x,a11))', ionPotArr2[0:2] )
     


; output strings with ionization potential to files

      printf, lunOut1, outStr1
      printf, lunOut2, outStr2

   endif
   
endwhile


; close input/output files

free_lun, lunOut2
free_lun, lunOut1

free_lun, lunIn

end
