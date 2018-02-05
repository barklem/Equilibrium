#!/bin/csh

set data_dir='/Users/remo/Desktop/barklem/data'
set levels_dir=${data_dir}'/levels_nist'
set current_dir=$PWD

cd ${levels_dir}

foreach elem ( 'H' 'He' 'Li'  'Be'  'B'   'C'   'N'  'O'  'F'  'Ne'   \
               'Na'  'Mg'  'Al'  'Si'  'P'  'S'  'Cl'  'Ar'  \
	       'K'   'Ca'  'Sc'  'Ti'  'V'  'Cr'  'Mn'  'Fe' \
               'Co'  'Ni'  'Cu'  'Zn'  'Ga'  'Ge'  'As'  'Se'  'Br'  'Kr'  \
	       'Rb'  'Sr'  'Y'   'Zr'  'Nb'  'Mo'  'Tc'  'Ru' \
               'Rh'  'Pd'  'Ag'  'Cd'  'In'  'Sn'  'Sb'  'Te'  'I'  'Xe'  \
	       'Cs'  'Ba'  'La'  'Ce'  'Pr'  'Nd'  'Pm'  'Sm'  \
               'Eu'  'Gd'  'Tb'  'Dy'  'Ho'  'Er'  'Tm'  'Yb'  \
               'Lu'  'Hf'  'Ta'  'W'  'Re'  'Os'  'Ir'  'Pt'  \
               'Au'  'Hg'  'Tl' 'Pb'  'Bi'  'Po'  'At'  'Rn' \
	       'Fr'  'Ra'  'Ac'  'Th'  'Pa'  'U' )
	       
	foreach ion ( 'I' 'II' 'III' )	
	        echo ${elem} ${ion}
		wget "http://physics.nist.gov/cgi-bin/ASD/energy1.pl?spectrum="${elem}"+"${ion}"&units=1&level_out=checked&j_out=checked&format=1" -O nist.tmp.0
		cp nist.tmp.0 ${elem}'_'${ion}'_nist.html'
		sed -e :a -e 's/<[^>]*>//g;/</N;//ba'  nist.tmp.0 > nist.tmp.1
		grep "|" nist.tmp.1 > nist.tmp.2
		sed -e '/Level/s/|//g'  nist.tmp.2 > nist.tmp.3
		sed -e 's/Level/E/g' nist.tmp.3 > ${elem}'_'${ion}'_nist.atom'
	end

end

cd ${current_dir}
