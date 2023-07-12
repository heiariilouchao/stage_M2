CC = gcc
OPT = -O3 -Wall
BIN_DIR = bin
BUILD_DIR = build


rdf: parse.o read.o rdf.o
	$(CC) $(OPT) -o $(BIN_DIR)/rdf $(addprefix $(BUILD_DIR)/,$?) -lm

msd: parse.o read.o msd.o
	$(CC) $(OPT) -o $(BIN_DIR)/msd $(addprefix $(BUILD_DIR)/,$?) -lm -fopenmp

# %.o: %/%.c
# 	$(CC) $(OPT) -o $(BUILD_DIR)/$@ -c $<

# Is equivalent to all these
parse.o: parse/parse.c
	$(CC) $(OPT) -o $(BUILD_DIR)/parse.o -c parse/parse.c

read.o: read/read.c
	$(CC) $(OPT) -o $(BUILD_DIR)/read.o -c read/read.c

rdf.o: rdf/rdf.c
	$(CC) $(OPT) -o $(BUILD_DIR)/rdf.o -c rdf/rdf.c

msd.o: msd/msd.c
	$(CC) $(OPT) -o $(BUILD_DIR)/msd.o -c msd/msd.c -fopenmp
	

reset: clean
	rm $(BIN_DIR)/*

.PHONY: clean rebuild

clean:
	rm $(BUILD_DIR)/*.o

rebuild: reset build