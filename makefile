# make file for lwacm

CC := gcc
C_FLAG := -Wall -Wno-unused-variable #-Wno-unused-label

all: lwacm

lwacm: lwacm.c
	$(CC) lwacm.c -o lwacm  $(C_FLAG)
	
testf: testf.c
	$(CC) testf.c -o testf  $(C_FLAG)

clean:
	rm -rf lwacm lwacm.c~ makefile~ README.md~ testf testf.c~ 
