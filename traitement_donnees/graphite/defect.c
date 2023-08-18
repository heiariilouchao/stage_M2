#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <errno.h>

#include "../utils/utils.h"
#include "../parse/parse.h"
#include "../read/read.h"
#include "../rdf/rdf.h"
#include "carbons.h"

char doc[] = "Processing and extracting the graphite informations related to the defects.";

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

	/* Extracting the carbons */
	int *N_carbons;
	Carbon **carbons;
	if ((errno = extract_carbons(N_configurations, N_atoms, atoms, box, &N_carbons, &carbons)) != 0)
		goto READ;

	/* Computing the layers */
	const double sep = 3.5;
	compute_layers(N_configurations, box, N_carbons, &carbons, sep);

	/* Computing the charges histogram */
	// The electrodes parameter: indicating the z-coordinate of the lower electrode
	const double limits[2] = {box[0].z_min, box[0].z_min + (box[0].z_max - box[0].z_min) / 2.};
	Electrode electrodes[2] = {LowerElectrode, UpperElectrode};

	// Computing the electrodes
	compute_electrodes(N_configurations, N_carbons, &carbons, limits, electrodes);

    // Converting the Carbons into Atoms
    Atom **carbons_atoms;
    if ((errno = convert_carbons(N_configurations, N_carbons, carbons, &carbons_atoms)) != 0)
        goto EXTRACT_CARBONS;

    /* Computing the RDFs */
    // Selecting the direct neighbors of the defect
    int *N_sp;
    Atom **sp;
    if ((errno = select_valency(N_configurations, N_carbons, LowerEqual, 2, carbons_atoms, &N_sp, &sp)) != 0)
        goto CONVERT_CARBONS;
    
    // Selecting the sodium ions
    int *N_sodium;
    Atom **sodium;
    if ((errno = select_elements(N_configurations, N_atoms, "Na", atoms, &N_sodium, &sodium)) != 0)
        goto SELECT_SP;
    
    // Computing the RDF between the two groups
    double *r, *rdf;
    if ((errno = compute_rdf(N_configurations, box, 100, false, N_sp, sp, N_sodium, sodium, &r, &rdf)) != 0)
        goto SELECT_SODIUM;
    
    if ((errno = write_rdf("output/graphite-defect/rdf_sp-Na.dat", "sp,Na", 100, r, rdf)) != 0)
        goto COMPUTE_RDF;
    
    // Computing the charges near the defect
    Group group;
    if ((errno = compute_average(N_configurations, N_sp, sp, Q, &group, "Average charge of the defect's neighbors")) != 0)
        goto COMPUTE_RDF;
    
    if ((errno = write_average("output/graphite-defect/q_defect.dat", N_configurations, steps, group)) != 0)
        goto COMPUTE_AVERAGE;

    /* Success */
    exit(EXIT_SUCCESS);

    /* Errors */
    COMPUTE_AVERAGE:
        free(group.N), free(group.average);
    COMPUTE_RDF:
        free(rdf), free(r);
    SELECT_SODIUM:
        for (int c = 0; c < N_configurations; c++)
            free(sodium[c]);
        free(sodium), free(N_sodium);
    SELECT_SP:
        for (int c = 0; c < N_configurations; c++)
        {
            for (int a = 0 ; a < N_sp[c] ; a++)
                free(sp[c][a].bonded);
            
            free(sp[c]);
        }
        free(sp), free(N_sp);
    CONVERT_CARBONS:
        for (int c = 0 ; c < N_configurations ; c++)
        {
            for (int a = 0 ; a < N_carbons[c] ; a++)
                free(carbons_atoms[c][a].bonded);
            
            free(carbons_atoms[c]);
        }
        free(carbons_atoms);
    EXTRACT_CARBONS:
        for (int c = 0; c < N_configurations; c++)
        {
            for (int a = 0 ; a < N_carbons[c] ; a++)
                free(carbons[c][a].atom.bonded);

            free(carbons[c]);
        }
        free(carbons), free(N_carbons);
    READ:
        free(box);
        for (int c = 0; c < N_configurations; c++)
            free(atoms[c]);
        free(atoms), free(N_atoms), free(steps);
    ERROR:
        exit(EXIT_FAILURE);
}