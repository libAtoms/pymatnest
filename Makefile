RngStream.so: RngStream.c RngStream.h
	gcc -O -fPIC -shared RngStream.c -o RngStream.so

MC_MD_steps.pdf: MC_MD_steps.tex ns_run
	egrep '#DOC' ns_run | sed 's/#DOC[ ]*//' > MC_MD_steps_body.tex
	pdflatex MC_MD_steps.tex
