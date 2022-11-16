#!/bin/csh

set data_dir='/Users/remo/Desktop/barklem/data'
set levels_dir=${data_dir}'/levels_nist'
set current_dir=$PWD

cd ${levels_dir}

foreach elem ( 'H' 'D' 'He')
	       
	foreach ion ( 'I' 'II' 'III' )	
	        echo ${elem} ${ion}
		wget "http://physics.nist.gov/cgi-bin/ASD/energy1.pl?spectrum="${elem}"+"${ion}"&level_out=checked&j_out=checked&format=1" -O nist.tmp.0
		cp nist.tmp.0 ${elem}'_'${ion}'_nist.html'
		sed -e :a -e 's/<[^>]*>//g;/</N;//ba'  nist.tmp.0 > nist.tmp.1
		grep "|" nist.tmp.1 > nist.tmp.2
		sed -e '/Level/s/|//g'  nist.tmp.2 > nist.tmp.3
		sed -e 's/Level/E/g' nist.tmp.3 > ${elem}'_'${ion}'_nist.atom'
	end

end

cp ../empty_nist.atom H_II_nist.atom
cp ../empty_nist.atom D_II_nist.atom
cp ../empty_nist.atom H_III_nist.atom
cp ../empty_nist.atom D_III_nist.atom
cp ../empty_nist.atom He_III_nist.atom

cd ${current_dir}
