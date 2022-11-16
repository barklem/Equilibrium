# DATA

Data sets:

"CDS format data" contains the original published data in ascii formats put on CDS, plus versions of the same with the November 2022 minor bug fix (vNov2022).  This is the same as published on github at public-data/partition-functions_and_equilibrium-constants

Then I provide various data sets in idl binary format and/or ascii.  The main reason for this is to provide calculations done on wider and more dense temperature grids, especially at high T.  But this should not imply the data are good at these temperatures - it is mostly to aid interpolation and extrapolation in code, where for practical applications something has to be calculated irrespective of whether it is good physically or not.

I have used four grids with silly names:

; Temperature grid for Barklem & Collet (2016) paper ("published T grid")<br>

T = [1e-5,  1e-4,  1e-3,  1e-2,  0.1,   0.15,  0.2,   0.3,   0.5,   0.7, $<br>
     1.0,   1.3,   1.7,   2.0,   3.0,   5.0,   7.,    10.,   15.,   20., $<br>
     30.,   50.,   70.,   100.,  130.,  170.,  200.,  250.,  300.,  500., $<br>
     700.,  1000., 1500., 2000., 3000., 4000., 5000., 6000., 7000., 8000., $<br>
     9000., 10000. ]<br>

; Temperature grid, Oct. 2016 calculations  ("extended T grid")<br>
  T = 10.^(findgen(1001)/100.-6)<br>

; Temperature grid, Aug. 2020 calculations, 10^-6 to 10^6, 100 or 400 points per decade<br>
  T = 10.^(findgen(12*100+1)/100.-6)        ;   ("vlarge T grid")<br>
  T = 10.^(findgen(12*400+1)/400.-6)        ;   ("extremely large T grid")<br>

  Inside the folders you can always find idl binary files:<br>

  allatom_collated.idl     - atomic partition functions<br>
  allequil_collated.idl    - molecular equilibrium constants<br>
  allpartf_collated.idl    - molecular partition functions<br>

  In some cases I also give the ascii versions of these tables, respectively:<br>

  atompartf_table.txt<br>
  equil_table.txt<br>
  partf_table.txt<br>
