F90_MODELS = $(wildcard *_model.F90)

all: doc RngStream.so fortran_MC_MD.so $(F90_MODELS:.F90=.so)

doc: MC_MD_steps.pdf

MC_MD_steps.pdf: MC_MD_steps.tex ns_run
	egrep '#DOC' ns_run | sed 's/#DOC[ ]*//' > MC_MD_steps_body.tex
	pdflatex MC_MD_steps.tex

# extra dependency for RngStream
RngStream.so: RngStream.c RngStream.h

# rules for making .so from .F90 and from .c
%.so : %.F90
	gfortran -g -O3 -fPIC -shared $< -o $@

%.so : %.c
	gcc -g -O3 -fPIC -shared $< -o $@

-include Makefile.model
