make clean
make html
make latex
pushd build/latex
# Uncomment when latexmk available
#make all-pdf
latex miri_tools.tex
latex miri_tools.tex
pdflatex miri_tools.tex
popd
