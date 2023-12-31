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
#include "carbons.h"

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
	Carbon **carbons;
	if ((errno = extract_carbons(N_configurations, N_atoms, atoms, box, &N_carbons, &carbons)) != 0)
		goto READ;
	
	// Computing the layers
	const double sep = 3.5;
	compute_layers(N_configurations, box, N_carbons, &carbons, sep);
	
	// Selecting the outer layers
	int *N_outer;
	Carbon **outer;
	if ((errno = select_layer(N_configurations, N_carbons, carbons, Outer, &N_outer, &outer)) != 0)
		goto EXTRACT_CARBONS;
	
	// Converting the outer Carbons to Atoms
	Atom **outer_atoms;
	if ((errno = convert_carbons(N_configurations, N_outer, outer, &outer_atoms)) != 0)
		goto SELECT_OUTER;

	// Selecting the lower electrode
	int *N_lower;
	Atom **lower;
	if ((errno = select_coordinate(N_configurations, N_outer, Lower, Coord_Z, box[0].z_min + (box[0].z_max - box[0].z_min) / 2., outer_atoms, &N_lower, &lower)) != 0)
		goto CONVERT_OUTER;


	/* Selecting the sodium ions */
	int *N_sodium;
	Atom **sodium;
	if ((errno = select_elements(N_configurations, N_atoms, "Na", atoms, &N_sodium, &sodium)) != 0)
		goto SELECT_LOWER;
	

	/* Selecting the upper carbons */
	int *N_upper;
	Atom **upper;
	if ((errno = select_coordinate(N_configurations, N_outer, Greater, Coord_Z, box[0].z_min + (box[0].z_max - box[0].z_min) / 2., outer_atoms, &N_upper, &upper)) != 0)
		goto RDF;
	

	/* Selecting the hydroxide ions */
	int *N_HO;
	Atom **HO;
	// Selecting the hydrogens and oxygens
	if ((errno = select_elements(N_configurations, N_atoms, "H,O", atoms, &N_HO, &HO)) != 0)
		goto SELECT_UPPER;

	// Computing the bonds
	if ((errno = compute_cutoff_bonds(N_configurations, N_HO, box, &HO, 1.1)) != 0)
		goto SELECT_HO;
	
	// Selecting the oxygens
	int *N_O;
	Atom **O;
	if ((errno = select_elements(N_configurations, N_HO, "O", HO, &N_O, &O)) != 0)
		goto SELECT_HO;
	
	// Selecting the hydroxide oxygens
	int *N_hydroxide;
	Atom **hydroxide;
	if ((errno = select_valency(N_configurations, N_O, LowerEqual, 1, O, &N_hydroxide, &hydroxide)) != 0)
		goto SELECT_O;
	

	/* Computing the RDFs */
	double *r, *rdf;

	// Between the lower carbons and the sodium ions
	if ((errno = compute_rdf(N_configurations, box, 100, false, N_lower, lower, N_sodium, sodium, &r, &rdf)) != 0)
		goto SELECT_SODIUM;

	if ((errno = write_rdf("output/graphite-rdf/rdf_lower-Na.dat", "lower,Na", 100, r, rdf)) != 0)
		goto RDF;

	free(r), free(rdf);

	// Between the lower carbons and the hydroxide ions
	if ((errno = compute_rdf(N_configurations, box, 100, false, N_lower, lower, N_hydroxide, hydroxide, &r, &rdf)) != 0)
		goto SELECT_SODIUM;

	if ((errno = write_rdf("output/graphite-rdf/rdf_lower-OH.dat", "lower,OH", 100, r, rdf)) != 0)
		goto RDF;

	free(r), free(rdf);

	// Between the upper carbons and the sodium ions
	if ((errno = compute_rdf(N_configurations, box, 100, false, N_upper, upper, N_sodium, sodium, &r, &rdf)) != 0)
		goto SELECT_HYDROXIDE;
	
	if ((errno = write_rdf("output/graphite-rdf/rdf_upper-Na.dat", "upper,Na", 100, r, rdf)) != 0)
		goto SELECT_HYDROXIDE;
	
	free(r), free(rdf);

	// Between the upper carbons and the hydroxide ions
	if ((errno = compute_rdf(N_configurations, box, 100, false, N_upper, upper, N_hydroxide, hydroxide, &r, &rdf)) != 0)
		goto SELECT_HYDROXIDE;
	
	if ((errno = write_rdf("output/graphite-rdf/rdf_upper-OH.dat", "upper,OH", 100, r, rdf)) != 0)
		goto SELECT_HYDROXIDE;
    

    /* Exiting normally */
	exit(EXIT_SUCCESS);


    /* Error handling */
	SELECT_HYDROXIDE:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0; a < N_hydroxide[c]; a++)
				free(hydroxide[c][a].bonded);
			
			free(hydroxide[c]);
		}
		free(hydroxide), free(N_hydroxide);
	SELECT_O:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0; a < N_O[c]; a++)
				free(O[c][a].bonded);
			
			free(O[c]);
		}
		free(O), free(N_O);
	SELECT_HO:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0; a < N_HO[c]; a++)
				free(HO[c][a].bonded);
			
			free(HO[c]);
		}
		free(HO), free(N_HO);
	SELECT_UPPER:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0; a < N_upper[c]; a++)
				free(upper[c][a].bonded);
			
			free(upper[c]);
		}
		free(upper), free(N_upper);
	RDF:
		free(r), free(rdf);
	SELECT_SODIUM:
		for (int c = 0 ; c < N_configurations ; c++)
			free(sodium[c]);
		free(sodium), free(N_sodium);
	SELECT_LOWER:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0; a < N_lower[c]; a++)
				free(lower[c][a].bonded);
			
			free(lower[c]);
		}
		free(lower), free(N_lower);
	
	CONVERT_OUTER:
		for (int c = 0; c < N_configurations; c++)
		{
			for (int a = 0; a < N_outer[c] ; a++)
				free(outer_atoms[c][a].bonded);
			
			free(outer_atoms[c]);
		}
		free(outer_atoms);
	SELECT_OUTER:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0 ; a < N_outer[c] ; a++)
				free(outer[c][a].atom.bonded);
			
			free(outer[c]);
		}
		free(outer), free(N_outer);
	EXTRACT_CARBONS:
		for (int c = 0 ; c < N_configurations ; c++)
		{
			for (int a = 0 ; a < N_carbons[c] ; a++)
				free(carbons[c][a].atom.bonded);
			
			free(carbons[c]);
		}
		free(carbons), free(N_carbons);
	READ:
		for (int c = 0 ; c < N_configurations ; c++)
			free(atoms[c]);
		free(atoms);
		free(box), free(N_atoms), free(steps);
	ERROR: exit(EXIT_FAILURE);
}
