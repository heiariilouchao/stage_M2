CC = gcc
OPT = -O3 -Wall
BIN_DIR = bin
BUILD_DIR = build


build: bonds.o parse.o process_graphite.o read.o
	$(CC) $(OPT) -o $(BIN_DIR)/process_graphite $(BUILD_DIR)/*.o

%.o: %.c
	$(CC) $(OPT) -o $(BUILD_DIR)/$@ -c $<

reset: clean
	rm $(BIN_DIR)/process_graphite

.PHONY: clean rebuild

clean:
	rm $(BUILD_DIR)/*.o

rebuild: reset build
