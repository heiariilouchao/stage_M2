# ifndef READ_H
# define READ_H


# include <stdio.h>

# include "utils.h"


# define STR_BUFF_LIMIT 256
# define INITIAL_CONF 10
# define INCREMENT_CONF 10


int read_trajectory(char *file_name, int timestep, int N_elements, char *labels, char **elements, int *N_conf, int **steps, int **N_selection, double ***bounds, Atom ***atoms);

int select_atom(FILE *file, char *labels, char *element);


# endif