CC = gcc
OPT = -g -O3 -Wall -Werror
BIN_DIR = bin
BUILD_DIR = build


rdf: utils.o parse.o read.o rdf.o
	$(CC) $(OPT) -o $(BIN_DIR)/rdf $(addprefix $(BUILD_DIR)/,$?) -lm

msd: utils.o parse.o read.o msd.o
	$(CC) $(OPT) -o $(BIN_DIR)/msd $(addprefix $(BUILD_DIR)/,$?) -lm -fopenmp

graphite: utils.o parse.o read.o bonds.o carbons.o graphite.o
	$(CC) $(OPT) -o $(BIN_DIR)/graphite $(addprefix $(BUILD_DIR)/,$?) -lm

graphite-rdf: utils.o parse.o read.o bonds.o rdf.o carbons.o graphite-rdf.o
	$(CC) $(OPT) -o $(BIN_DIR)/graphite-rdf $(addprefix $(BUILD_DIR)/,$?) -lm

defect: utils.o parse.o read.o rdf.o bonds.o carbons.o defect.o
	$(CC) $(OPT) -o $(BIN_DIR)/defect $(addprefix $(BUILD_DIR)/,$?) -lm

# Is equivalent to all these
utils.o: utils/utils.c
	$(CC) $(OPT) -o $(BUILD_DIR)/utils.o -c utils/utils.c

parse.o: parse/parse.c
	$(CC) $(OPT) -o $(BUILD_DIR)/parse.o -c parse/parse.c

read.o: read/read.c
	$(CC) $(OPT) -o $(BUILD_DIR)/read.o -c read/read.c

rdf.o: rdf/rdf.c
	$(CC) $(OPT) -o $(BUILD_DIR)/rdf.o -c rdf/rdf.c

msd.o: msd/msd.c
	$(CC) $(OPT) -o $(BUILD_DIR)/msd.o -c msd/msd.c -fopenmp

bonds.o: bonds/bonds.c
	$(CC) $(OPT) -o $(BUILD_DIR)/bonds.o -c bonds/bonds.c

carbons.o: graphite/carbons.c
	$(CC) $(OPT) -o $(BUILD_DIR)/carbons.o -c graphite/carbons.c

graphite.o: graphite/graphite.c
	$(CC) $(OPT) -o $(BUILD_DIR)/graphite.o -c graphite/graphite.c

graphite-rdf.o: graphite/rdf.c
	$(CC) $(OPT) -o $(BUILD_DIR)/graphite-rdf.o -c graphite/rdf.c

defect.o: graphite/defect.c
	$(CC) $(OPT) -o $(BUILD_DIR)/defect.o -c graphite/defect.c

reset: clean
	rm $(BIN_DIR)/*

.PHONY: clean rebuild

clean:
	rm $(BUILD_DIR)/*.o
