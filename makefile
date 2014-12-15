# make file for lwacm

CC := gcc
C_FLAG := -Wall -Wno-unused-variable #-Wno-unused-label

all: lwacm

.PHONEY: test

lwacm: lwacm.c
	$(CC) lwacm.c -o lwacm  $(C_FLAG)
	
test: test.c
	$(CC) test.c -o test  $(C_FLAG)

clean:
	rm -rf lwacm lwacm.c~ makefile~ README.md~ test test.c~ 
