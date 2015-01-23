# make file for lwacm
.PHONY: test

CC := gcc
CFLAGS := -Wall -O3

SRC := $(wildcard ./*.c)
BIN := $(SRC:.c=)


default: $(BIN)

test: mem_test

clean:
	rm -rf $(BIN) *.c~ lwacm_raw~ makefile~ README.md~ lwacm_run_* mem_test
	
	
