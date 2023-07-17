include Makefile.arch

F90_MODEL_UTILS = example_mat_mod.F90 example_LJ_params_mod.F90 example_bead_spring_polymer_params_mod.F90 example_Jagla_params_mod.F90 example_CollapsingSpheres_params_mod.F90 Wang_LJ_params_mod.F90

F90_MODELS = $(wildcard *_model.F90)

all: RngStream.so libs

libs: fortran_MC_MD.so $(F90_MODELS:.F90=.so)

#doc: MC_MD_steps.pdf

#MC_MD_steps.pdf: MC_MD_steps.tex ns_run.py
#	egrep '#DOC' ns_run.py | sed 's/#DOC[ ]*//' > MC_MD_steps_body.tex
#	pdflatex MC_MD_steps.tex

# extra dependency for RngStream
RngStream.so: RngStream.c RngStream.h

# rules for making .o and .so from .F90 and from .c
%.so : %.F90 $(F90_MODEL_UTILS:.F90=.o) 
	$(FC) $(FFLAGS) $(SO_FLAGS) $< $(F90_MODEL_UTILS:.F90=.o) -o $@

%.o : %.F90
	$(FC) $(FFLAGS) -c $< -o $@

%.so : %.c
	$(CC) $(CFLAGS) $(SO_FLAGS) $< -o $@

-include Makefile.model
