# ifndef UTILS_H
# define UTILS_H


# define STR_GROUP_DESC_LIMIT 256

typedef struct Atom
{
    int serial;
    double x, y, z;
    double q;

    int N_bonds;
    int *bonded;
    
    int group;
} Atom;


typedef struct Group
{
    int N;
    double average;
} Group;

# endif