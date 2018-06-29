make clean
make html
make latex
pushd build/latex
#make all-pdf
latex miri_simulators.tex
latex miri_simulators.tex
pdflatex miri_simulators.tex
popd
