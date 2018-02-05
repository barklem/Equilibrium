PRO default, var, val, set=set
   COMMON cdefault, res
   if n_elements(set) ne 0 then res=set
   if n_elements(val) eq 0 then return
   if n_elements(res) eq 0 then res=0
   if n_elements(var) eq 0 or res ne 0 then var=val
END