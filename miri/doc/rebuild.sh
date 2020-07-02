make clean
make html
make latex
pushd build/latex
#make all-pdf
latex miridoc.tex
latex miridoc.tex
pdflatex miridoc.tex
popd
