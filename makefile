# make file for lwacm

CC := gcc
CFLAGS := -Wall -O3 -static

SRC := $(wildcard ./*.c)
BIN := $(SRC:.c=)


default: $(BIN)

clean:
	rm -rf $(BIN)


