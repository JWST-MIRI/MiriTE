make clean
make html
make latex
pushd build/latex
# Uncomment when latexmk available
#make all-pdf
latex miri_datamodels.tex
latex miri_datamodels.tex
pdflatex miri_datamodels.tex
popd
