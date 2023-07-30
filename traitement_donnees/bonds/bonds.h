# ifndef BONDS_H
# define BONDS_H

# include "../utils/utils.h"

int compute_cutoff_bonds(int N_conf, int *N_selection, Box *box, Atom ***atoms, const double R);

# endif