# ifndef BONDS_H
# define BONDS_H

# include "../utils.h"

int compute_cutoff_bonds(int N_conf, int *N_selection, double **bounds, Atom ***atoms, const double R);

# endif