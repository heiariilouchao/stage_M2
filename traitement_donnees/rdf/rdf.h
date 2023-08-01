# ifndef RDF_H
# define RDF_H


# include "../utils/utils.h"


int compute_pairs(int N_elements, char **elements, int *N_pairs, char ***pairs);


double compute_cutoff(Box *box);


int compute_rdf(int N_configurations, Box *box, int N_bins, bool are_identical, int *N1, Atom **a1, int *N2, Atom **a2, double **r, double **rdf);


int write_pairs(char *file_name, int N_pairs, char **pairs, int N_bins, double *r, double **rdf);


int write_rdf(char *file_name, char *label, int N_bins, double *r, double *rdf);

# endif