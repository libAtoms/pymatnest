all: doc fortran_MC_MD.so example_model.so

doc: MC_MD_steps.pdf

fortran_MC_MD.so: fortran_MC_MD.F90
	gfortran -g -O -fPIC -shared fortran_MC_MD.F90 -o fortran_MC_MD.so

example_model.so: example_model.F90
	gfortran -g -O -fPIC -shared example_model.F90 -o example_model.so

MC_MD_steps.pdf: MC_MD_steps.tex ns_run
	egrep '#DOC' ns_run | sed 's/#DOC[ ]*//' > MC_MD_steps_body.tex
	pdflatex MC_MD_steps.tex

RngStream.so: RngStream.c RngStream.h
	gcc -O -fPIC -shared RngStream.c -o RngStream.so

-include Makefile.model
