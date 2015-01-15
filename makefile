# make file for lwacm

CC := gcc
C_FLAG := -Wall -Wno-unused-variable #-Wno-unused-label

all: lwacm  lwacm.test

lwacm: lwacm.c
	$(CC) lwacm.c -o lwacm  $(C_FLAG)
	
lwacm.test: lwacm.test.c
	$(CC) lwacm.test.c -o lwacm.test  $(C_FLAG)

clean:
	rm -rf lwacm lwacm.c~ makefile~ README.md~ lwacm.test lwacm.test.c~ 
