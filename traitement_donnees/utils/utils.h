# ifndef UTILS_H
# define UTILS_H


# define STR_ELEMENT_LIMIT 3
# define STR_GROUP_DESCRIPTION_LIMIT 256
# define MAXIMUM_N_BONDS 4


typedef struct Box
{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
} Box;


typedef struct Atom
{
    int serial;
    char element[STR_ELEMENT_LIMIT];

    double x, y, z;
    double xu, yu, zu;
    double q;

    int N_bonds;
    int *bonded;

    int group;
} Atom;


typedef struct Group
{
    char description[STR_GROUP_DESCRIPTION_LIMIT];

    int *N;
    double *average;
} Group;


typedef enum ComparisonOperator
{
    Greater,
    GreaterEqual,
    Equal,
    LowerEqual,
    Lower
} ComparisonOperator;


typedef enum Coordinate
{
    X,
    Y,
    Z
} Coordinate;


int select_elements(int N_configurations, int *N_atoms, char *labels, Atom **all, int **N_selected, Atom ***selected);


int select_valency(int N_configurations, int *N_atoms, ComparisonOperator operator, int valency, Atom **all, int **N_selected, Atom ***selected);


int select_coordinate(int N_configurations, int *N_atoms, ComparisonOperator operator, Coordinate coordinate, double value, Atom **all, int **N_selected, Atom ***selected);


int average_charge(int N_configurations, int *N_atoms, Atom **atoms, Group *group, char *description);


int write_average(char *file_name, int N_configurations, int *steps, Group group);

# endif