# make file for lwacm

CC := gcc
C_FLAG := -Wall -Wno-unused-variable #-Wno-unused-label

all: lwacm file_gen

lwacm: lwacm.c
	$(CC) lwacm.c -o lwacm  $(C_FLAG)
file_gen: file_gen.c
	$(CC) file_gen.c -o file_gen  $(C_FLAG)
clean:
	rm -rf lwacm lwacm.c~ makefile~ README.md~ lwacm.test
