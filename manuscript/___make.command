#!/bin/bash


echo -e "\n\n~~~~~~~~~~~~~~~~~~\n"
echo Creating the LaTeX document...
echo -e "\n~~~~~~~~~~~~~~~~~~\n"

## This creates the PDF from the __ms.tex document.
## On macOS, you should just be able to click on it from Finder

export name="charr_ms"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR"

pdflatex -interaction=nonstopmode ${name}
bibtex ${name}
pdflatex -interaction=nonstopmode ${name}
pdflatex -interaction=nonstopmode ${name}
rm -Rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.synctex.gz



exit