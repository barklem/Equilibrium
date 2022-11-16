#!/bin/csh

# The script downloads tables with energy level information for some elements (first and second ionization stages)

set data_dir='/Users/remo/Desktop/barklem/data'
set tables_dir=${data_dir}'/tables_nist'
set current_dir=$PWD

cd ${tables_dir}

foreach element ( 'astatine' 'protactinium' 'uranium' 'thorium' 'radon' 'polonium' 'francium' )
	       
    wget "http://physics.nist.gov/PhysRefData/Handbook/Tables/"${element}"table5.htm" -O table5.html
    wget "http://physics.nist.gov/PhysRefData/Handbook/Tables/"${element}"table6.htm" -O table6.html
    grep "table7.htm" table5.html > ${element}"_I_table".html
    grep "table7.htm" table6.html > ${element}"_II_table".html
# sed -e :a -e 's/<[^>]*>//g;/</N;//ba'  nist.tmp.0 > nist.tmp.1
#		grep "|" nist.tmp.1 > nist.tmp.2
#		sed -e '/Level/s/|//g'  nist.tmp.2 > nist.tmp.3
#		sed -e 's/Level/E/g' nist.tmp.3 > ${elem}'_'${ion}'_nist.atom'

end

cd ${current_dir}
