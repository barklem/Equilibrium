#!/bin/csh
foreach i (*_*_nist.atom)
  sed -e '/J/d' $i > tmp
  cat tmp endline > tmp2
  mv tmp2 $i
end
