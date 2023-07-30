# ifndef READ_H
# define READ_H


# include <stdio.h>

# include "../utils/utils.h"


# define STR_BUFF_LIMIT 256
# define INITIAL_CONF 10
# define INCREMENT_CONF 10


int read_elements(Arguments **args);

int read_trajectory(Arguments *arguments, int *N_conf, int **steps, int **N_selection, Box **box, Atom ***atoms);

int select_atom(FILE *file, char *labels, char *element);


# endif