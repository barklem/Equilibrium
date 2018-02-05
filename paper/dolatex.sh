#!/bin/bash
export texfile='equil'
pdflatex ${texfile}
bibtex ${texfile}
pdflatex ${texfile}
pdflatex ${texfile}
