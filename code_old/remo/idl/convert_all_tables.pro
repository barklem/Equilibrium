elem_names = [ 'astatine' , 'protactinium', 'uranium', 'thorium', 'radon', 'polonium', 'francium' ]
elem_syms  = [ 'At', 'Pa', 'U', 'Th', 'Rn', 'Po', 'Fr' ]

nelem = n_elements( elem_names )

for i=0L,nelem-1L do begin

   fin  = elem_names[i]+'_I_table.html'
   fout = elem_syms[i]+'_I_nist.atom'
   convert_nist_table, fin, fout

   fin  = elem_names[i]+'_II_table.html'
   fout = elem_syms[i]+'_II_nist.atom'
   convert_nist_table, fin, fout

endfor

end
