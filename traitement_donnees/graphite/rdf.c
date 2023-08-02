# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <argp.h>

# include "../utils/utils.h"
# include "../parse/parse.h"
# include "../read/read.h"
# include "../bonds/bonds.h"
# include "../rdf/rdf.h"

char doc[] = "Processing and extracting the graphite RDF-related informations.";

char args_doc[] = "CONF_FILE";


int main(int argc, char **argv)
{
    /* Parsing the arguments */
	Arguments arguments;

	// Options' default values
	set_default_options(&arguments);

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
	{
		perror("Parsing");
		exit(EXIT_FAILURE);
	}

    /* Reading the carbons */
	// Declaring the arrays
	int N_configurations, *steps, *N_atoms;
	Atom **atoms;
	Box *box;

	// Reading the file
	if ((errno = read_trajectory(&arguments, &N_configurations, &steps, &N_atoms, &box, &atoms)) != 0)
		goto ERROR;

	
	/* Selecting the atoms */
	int *N_carbons;
	Atom **carbons;
	if ((errno = select_elements(N_configurations, N_atoms, "C", atoms, &N_carbons, &carbons)) != 0)
		goto READ;
	

	/* Selecting the upper electrode */
	int *N_upper;
	Atom **upper;
	if ((errno = select_coordinate(N_configurations, N_carbons, Greater, Z, (box[0].z_max - box[0].z_min) / 2., carbons, &N_upper, &upper)) != 0)
		goto SELECT_CARBONS;



    /* Computing the bonds */
	if ((errno = compute_cutoff_bonds(N_configurations, N_carbons, box, &carbons, 1.7)) != 0)
		goto SELECT_CARBONS;
    

    /* Separating the electrodes */

    

    /* Exiting normally */
	exit(EXIT_SUCCESS);


    /* Error handling */
	RDF:
		free(r), free(rdf);
	SELECT_SP:
		for (int c = 0 ; c < N_configurations ; c++)
			free(sp_carbons[c]);
		free(sp_carbons), free(N_sp);
	SELECT_HYDROXYDE:
		for (int c = 0 ; c < N_configurations ; c++)
			free(hydroxyde[c]);
		free(hydroxyde), free(N_hydroxyde);
	SELECT_OH1:
		for (int c = 0 ; c < N_configurations ; c++)
			free(OH1[c]);
		free(OH1), free(N_OH1);
	SELECT_OH:
		for (int c = 0 ; c < N_configurations ; c++)
			free(OH[c]);
		free(OH), free(N_OH);
	SELECT_CARBONS:
		for (int c = 0 ; c < N_configurations ; c++)
			free(carbons[c]);
		free(carbons), free(N_carbons);
	READ:
		for (int c = 0 ; c < N_configurations ; c++)
			free(atoms[c]);
		free(atoms);
		free(box), free(N_atoms), free(steps);
	ERROR: exit(EXIT_FAILURE);
}
