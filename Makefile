RngStream.so: RngStream.c RngStream.h
	gcc -O -fPIC -shared RngStream.c -o RngStream.so
