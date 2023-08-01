# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <errno.h>

# include "../parse/parse.h"
# include "../utils/utils.h"
# include "../read/read.h"
# include "rdf.h"


char doc[] = "Computing the RDFs of a .lammpstrj configurations file.";

char args_doc[] = "CONF_FILE";


int main(int argc, char **argv)
{
	/* Parsing the arguments */
	Arguments arguments;

	// Options' default values
	set_default_options(&arguments);

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
		goto EXIT;


	/* Reading the configurations */
	int N_configurations, *N_atoms, *steps;
	Box *box;
	Atom **all;

	// Reading the file
	if ((errno = read_trajectory(&arguments, &N_configurations, &steps, &N_atoms, &box, &all)) != 0)
		goto EXIT;
	

	/* Computing the pairs*/
	int N_pairs;
	char **pairs;
	if ((errno = compute_pairs(arguments.N_elements, arguments.elements, &N_pairs, &pairs)) != 0)
		goto READ;
	

	/* Computing the RDFs */
	// The parameters
	int N_bins = 100;
	double *r, **rdf;

	// Allocating the arrays
	if ((rdf = malloc(N_pairs * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (rdf)");
		goto PAIRS;
	}

	// Computing each RDF
	for (int p = 0 ; p < N_pairs ; p++)
	{
		int *N1, *N2;
		Atom **a1, **a2;
		char* pair = strdup(pairs[p]);

		char *e1 = strtok(pair, ",");
		if ((errno = select_elements(N_configurations, N_atoms, e1, all, &N1, &a1)) != 0)
		{
			perror("Selecting atoms (a1)");
			goto RDF;
		}

		char *e2 = strtok(NULL, ",");
		if ((errno = select_elements(N_configurations, N_atoms, e2, all, &N2, &a2)) != 0)
		{
			perror("Selecting atoms (a2)");
			goto RDF;
		}

		bool are_identical = (strcmp(e1, e2) == 0);

		free(r);
		if ((errno = compute_rdf(N_configurations, box, N_bins, are_identical, N1, a1, N2, a2, &r, &(rdf[p]))) != 0)
			goto RDF;
		
		for (int c = 0 ; c < N_configurations ; c++)
			free(a1[c]), free(a2[c]);
		free(a1), free(a2);
		free(N1), free(N2);
	}
	

	/* Writing the output */
	if ((errno = write_pairs("output/rdf.dat", N_pairs, pairs, N_bins, r, rdf)) != 0)
		goto RDF;


	/* Success */
	exit(EXIT_SUCCESS);


	/* Error handling */
	RDF:
		for (int p = 0 ; p < N_pairs ; p++) free(rdf[p]);
		free(rdf);
		free(r);
	PAIRS:
		for (int p = 0 ; p < N_pairs ; p++) free(pairs[p]);
		free(pairs);
	READ:
		for (int c = 0 ; c < N_configurations ; c++)
			free(all[c]);
		free(all);
		free(box);
		free(steps);
		free(N_atoms);
	EXIT:
		exit(EXIT_FAILURE);
}