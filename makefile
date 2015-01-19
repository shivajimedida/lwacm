# make file for lwacm

CC := gcc
C_FLAG := -Wall -Wno-unused-variable #-Wno-unused-label

SRC := $(wildcard ./*.c)
BIN := $(SRC:.c=)


default: $(BIN)

clean:
	rm -rf lwacm code_gen *.c~ makefile~ README.md~ $(BIN)
