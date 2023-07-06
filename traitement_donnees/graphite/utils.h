#ifndef UTILS_H
#define UTILS_H

typedef struct Atom
{
    int serial;
    double x, y, z;
    double q;

    int N_bonds;
    int *bonded;
} Atom;

#endif