function relabel, s
; converts label string s to something more useful
; needs repstr from astrolib

snew = repstr(s, '&quot;', '"')
snew = repstr(snew, '&#931;', 'Sigma')
snew = repstr(snew, '&#963;', '-sigma')
snew = repstr(snew, '&#928;', 'Pi')
snew = repstr(snew, '&#960;', '-pi')
snew = repstr(snew, '&#916;', 'Delta')
snew = repstr(snew, '&#948;', '-delta')
snew = repstr(snew, '&#934;', 'Phi')
snew = repstr(snew, '&#966;', '-phi')
snew = repstr(snew, '&#937;', 'Omega')
snew = repstr(snew, '&gt;', '>')
snew = repstr(snew, '<sup>', '^')
snew = repstr(snew, '</sup>', ' ')
snew = repstr(snew, '<sub>', '_')
snew = repstr(snew, '</sub>', ' ')
snew = repstr(snew, ' )', ')')
snew = remove_links(snew)
snew = strtrim(snew)

snew = repstr(snew, 'B ^l Sigma^+', 'B ^1 Sigma^+') ; fixes a specific error in BeS state
snew = repstr(snew, 'L (^1 Sigma^+ , 1Pi)', 'L ^1 Sigma^+') ; fixes a specific error in HBr state
snew = repstr(snew, 'X ^1 Sigma_g ^  1s-sigma^2', 'X ^1 Sigma_g ^+  1s-sigma^2') ; fixes a specific error in H2 state

return, snew
end