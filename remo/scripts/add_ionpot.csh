#!/bin/csh

set data_dir='/Users/remo/Desktop/barklem/data'
set ionpot_dir=${data_dir}'/ionpot'
set levels_dir=${data_dir}'/levels'
set lev_ion_dir=${data_dir}'/levels_ionpot'

set current_dir=$PWD

set ionpot_file=${ionpot_dir}'/ionpot_H_He_metals_moog2007.txt'

set count   = `awk '{count++; print count;}'  $ionpot_file`
set atomnum = `awk '{print $1;}' $ionpot_file`
set el      = `awk '{print $2;}' $ionpot_file`
set chi1    = `awk '{print $3;}' $ionpot_file`
set chi2    = `awk '{print $4;}' $ionpot_file`
set chi3    = `awk '{print $5;}' $ionpot_file`

echo $count

cd $levels_dir


foreach i ($count)
   
  echo $i $atomnum[$i] $el[$i] $chi1[$i] $chi2[$i] $chi3[$i] 

  echo $chi1[$i] > tmp.1
  echo $chi2[$i] > tmp.2
  echo $chi3[$i] > tmp.3

  cat tmp.1 $el[$i]'_I_nist.atom' ${data_dir}"/end_nist.atom" > ion.tmp.1
  mv -f ion.tmp.1 ${lev_ion_dir}'/'${el[$i]}'_I_nist.atom'

  cat tmp.2 $el[$i]'_II_nist.atom' ${data_dir}"/end_nist.atom" > ion.tmp.2
  mv -f ion.tmp.2 ${lev_ion_dir}'/'${el[$i]}'_II_nist.atom'

  cat tmp.3 $el[$i]'_III_nist.atom' ${data_dir}"/end_nist.atom" > ion.tmp.3
  mv -f ion.tmp.3 ${lev_ion_dir}'/'${el[$i]}'_III_nist.atom'

end

cd ${current_dir}
