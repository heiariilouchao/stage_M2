# ifndef UTILS_H
# define UTILS_H


typedef struct Atom
{
    int serial;
    int element_ID;

    double x, y, z;
    double q;

    int N_bonds;
    int *bonded;

    int group;
} Atom;


# endif