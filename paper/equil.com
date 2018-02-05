latex equil
latex equil
bibtex equil
latex equil
dvips -t a4 -o equil.ps equil
ps2pdf equil.ps

