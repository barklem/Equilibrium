       program check
       
       real t
       double precision eqk, part, pion
       
       t= 3000.
       print *, 't = ', t
       call molcon('H2+', t, 2, 0.d0, 0.d0, eqk, part, pion)
       print *, 'H2+', eqk, part, pion
       call molcon2('H2+', t, 2, 0.d0, 0.d0, eqk, part, pion)
       print *, 'H2+', eqk, part, pion
       
       call molcon('H2', t, 2, 0.d0, 0.d0, eqk, part, pion)
       print *, 'H2', eqk, part, pion
       call molcon2('H2', t, 2, 0.d0, 0.d0, eqk, part, pion)
       print *, 'H2', eqk, part, pion
       
       call molcon('H2-', t, 2, 0.d0, 0.d0, eqk, part, pion)
       print *, 'H2-', eqk, part, pion
       call molcon2('H2-', t, 2, 0.d0, 0.d0, eqk, part, pion)
       print *, 'H2-', eqk, part, pion
       
       end
       
