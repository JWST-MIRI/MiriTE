make clean
make html
make latex
pushd build/latex
# Uncomment when latexmk available.
#make all-pdf
latex miri_scasim.tex
latex miri_scasim.tex
pdflatex miri_scasim.tex
popd
