CC = gcc
OPT = -O3 -Wall
BIN_DIR = bin
BUILD_DIR = build


build: parse.o read.o compute_rdf.o
	$(CC) $(OPT) -o $(BIN_DIR)/compute_rdf $(BUILD_DIR)/*.o -lm

%.o: %.c
	$(CC) $(OPT) -o $(BUILD_DIR)/$@ -c $<

reset: clean
	rm $(BIN_DIR)/compute_rdf

.PHONY: clean rebuild

clean:
	rm $(BUILD_DIR)/*.o

rebuild: reset build