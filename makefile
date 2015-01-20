# make file for lwacm

CC := gcc
CFLAGS := -Wall -O3

SRC := $(wildcard ./*.c)
BIN := $(SRC:.c=)


default: $(BIN)

clean:
	rm -rf lwacm code_gen *.c~ lwacm_raw~ makefile~ README.md~ lwacm_run_*
	
	
