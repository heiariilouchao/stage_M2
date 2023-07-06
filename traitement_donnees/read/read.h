#ifndef READ_H
#define READ_H

#define STR_BUFF_LIMIT 256
#define INITIAL_CONF 10
#define INCREMENT_CONF 10

typedef struct Atom
{
    int serial;
    double x, y, z;
    double q;

    int N_bonds;
    int *bonded;
} Atom;

int read_trajectory(char *file_name, int timestep, int *N_conf, int **steps, int **N_selection, double ***bounds, Atom ***atoms);

#endif