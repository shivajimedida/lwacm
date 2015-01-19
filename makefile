# make file for lwacm

CC := gcc
C_FLAG := -Wall -O3

SRC := $(wildcard ./*.c)
BIN := $(SRC:.c=)


default: $(BIN)

clean:
	rm -rf lwacm code_gen *.c~ makefile~ README.md~ lwacm_run_*
	
	
