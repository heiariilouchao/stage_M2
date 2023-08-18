/**
 * @file graphite.c
 * @author H. Lou Chao (heiariilouchao@gmail.com)
 * @brief Extracts the graphite-related informations from a .lammpstrj file.
 * @version 0.1
 * @date 2023-08-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <argp.h>

#include "../utils/utils.h"
#include "../parse/parse.h"
#include "../read/read.h"
#include "carbons.h"

char doc[] = "Processing and extracting the graphite informations.";

char args_doc[] = "CONF_FILE";

int write(char *file_name, int N_conf, int *N_selection, int *steps, Box *box, Atom **atoms, bool **layers, bool **electrodes)
{
	printf("Writing output...\n");

	/* Opening the output file */
	FILE *output = fopen(file_name, "w");
	if (output == NULL)
	{
		perror("Opening the output file");
		return EIO;
	}

	/* Writing */
	for (int c = 0; c < N_conf; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);

		fprintf(output, "ITEM: TIMESTEP\n%d\n", steps[c]);

		fprintf(output, "ITEM: NUMBER OF ATOMS\n%d\n", N_selection[c]);

		fprintf(output, "ITEM: BOX BOUNDS pp pp pp\n%lf %lf %lf %lf %lf %lf\n",
				box[c].x_min, box[c].x_max,
				box[c].y_min, box[c].y_max,
				box[c].z_min, box[c].z_max);

		fprintf(output, "ITEM: ATOMS id element x y z q nb id1 id2 id3 id4 in electrode\n");
		for (int a = 0; a < N_selection[c]; a++)
		{
			fprintf(output, "%d %s", atoms[c][a].serial, "C");
			fprintf(output, " %lf %lf %lf", atoms[c][a].x, atoms[c][a].y, atoms[c][a].z);
			fprintf(output, " %lf %d", atoms[c][a].q, atoms[c][a].N_bonds);
			for (int b = 0; b < atoms[c][a].N_bonds; b++)
				fprintf(output, " %d", atoms[c][atoms[c][a].bonded[b]].serial);
			for (int b = atoms[c][a].N_bonds; b < 4; b++)
				fprintf(output, " N/A");
			fprintf(output, " %d %d\n", layers[c][a], electrodes[c][a]);
		}
	}

	/* Exiting successfully */
	// Closing the file
	fclose(output);

	// Exiting
	return 0;
}

int average_carbons(int N_configurations, int *N_carbons, Carbon **carbons, AtomAttribute attribute, Group *group, char *description)
{
	/* Transforming the data into atoms */
	// Allocating the array
	Atom **atoms;
	if ((atoms = malloc(N_configurations * sizeof(Atom *))) == NULL)
	{
		perror("Allocating an array (atoms)");
		goto NOMEM;
	}

	// Transfering the data
	for (int c = 0; c < N_configurations; c++)
	{
		if ((atoms[c] = malloc(N_carbons[c] * sizeof(Atom))) == NULL)
		{
			perror("Allocating an array slot (atoms[])");
			goto ATOMS;
		}
		for (int a = 0; a < N_carbons[c]; a++)
			atoms[c][a] = carbons[c][a].atom;
	}

	/* Computing the average */
	// Computing the local group
	Group group_local;
	if (compute_average(N_configurations, N_carbons, atoms, attribute, &group_local, description) != 0)
	{
		perror("Computing the group");
		goto ATOMS_CONF;
	}

	// Allocating the arrays for the group
	if (((*group).N = malloc(N_configurations * sizeof(int))) == NULL)
	{
		perror("Allocating an array (group.N)");
		goto GROUP_LOCAL;
	}

	if (((*group).average = malloc(N_configurations * sizeof(double))) == NULL)
	{
		perror("Allocating an array (group.average)");
		goto GROUP_N;
	}

	// Transfering the data
	strncpy((*group).description, description, STR_GROUP_DESCRIPTION_LIMIT - 1);
	(*group).description[STR_GROUP_DESCRIPTION_LIMIT - 1] = '\0';
	for (int c = 0; c < N_configurations; c++)
	{
		(*group).N[c] = group_local.N[c];
		(*group).average[c] = group_local.average[c];
	}

	/* Success */
	// Freeing the temporary variables
	free(group_local.average), free(group_local.N);
	for (int c = 0; c < N_configurations; c++)
		free(atoms[c]);
	free(atoms);
	return 0;

	/* Errors */
GROUP_N:
	free((*group).N);
GROUP_LOCAL:
	free(group_local.average), free(group_local.N);
ATOMS_CONF:
	for (int c = 0; c < N_configurations; c++)
		free(atoms[c]);
ATOMS:
	free(atoms);
NOMEM:
	return ENOMEM;
}

int compute_density_histograms(int N_conf, int *N_selection, Atom **atoms, Box *box, double ***z, double ***densities, const int N_bins)
{
	printf("Computing the density histogram...\n");

	double delta = (box[0].z_max - box[0].z_min) / N_bins;

	/* Allocating the array */
	if ((*z = malloc(N_conf * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (z)");
		return ENOMEM;
	}

	if ((*densities = malloc(N_conf * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (densities)");
		goto Z;
	}

	for (int c = 0; c < N_conf; c++)
	{
		if (((*z)[c] = malloc(N_bins * sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (z[])");
			goto DENSITIES;
		}

		if (((*densities)[c] = calloc(N_bins, sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (densities[])");
			goto DENSITIES;
		}

		double z0 = box[c].z_min;
		for (int b = 0; b < N_bins; b++)
			(*z)[c][b] = z0 + delta * b + delta / 2.;
	}

	/* Computing the densities */
	for (int c = 0; c < N_conf; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0; a < N_selection[c]; a++)
		{
			int bin = (atoms[c][a].z - box[c].z_min) / delta;
			if (bin < N_bins)
			{
				// N[c]++;
				(*densities)[c][bin]++;
			}
		}

		for (int b = 0; b < N_bins; b++)
			// (*densities)[c][b] /= (bounds[c][1] - bounds[c][0]) * (bounds[c][3] - bounds[c][2]) * delta;
			(*densities)[c][b] /= N_selection[c];
	}

	/* Success */
	return 0;

/* Errors */
DENSITIES:
	free(densities);
Z:
	free(z);
	// N: free(N);
	return ENOMEM;
}

int compute_average_density(int N_configurations, int N_bins, double **z, double **densities, double **average_z, double **average_density)
{
	/* Allocating the arrays */
	if ((*average_z = calloc(N_bins, sizeof(double))) == NULL)
	{
		perror("Allocating an array (average_z)");
		goto NOMEM;
	}

	if ((*average_density = calloc(N_bins, sizeof(double))) == NULL)
	{
		perror("Allocating an array (average_density)");
		goto Z;
	}

	for (int c = 0 ; c < N_configurations ; c++)
	{
		for (int b = 0 ; b < N_bins ; b++)
		{
			(*average_z)[b] += z[c][b];
			(*average_density)[b] += densities[c][b];
		}
	}

	for (int b = 0 ; b < N_bins ; b++)
	{
		(*average_z)[b] /= N_configurations;
		(*average_density)[b] /= N_configurations;
	}

	/* Success */
	return 0;

	/* Errors */
	Z:
		free(*average_density);
	NOMEM:
		return ENOMEM;
}

int write_densities_hist(char *file_name, int N_conf, int N_bins, int *timestep, double **z, double **densities)
{
	printf("Writing the densities histograms...\n");

	/* Opening the file */
	FILE *output;

	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}

	/* Writing the output */
	// Header
	fprintf(output, "# step z density\n");

	// Data
	for (int c = 0; c < N_conf; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int b = 0; b < N_bins; b++)
			fprintf(output, "  %d %lf %lf\n", timestep[c], z[c][b], densities[c][b]);
		fprintf(output, "\n\n");
	}

	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}

int write_average_density(char *file_name, int N_bins, double *z, double *density)
{
	/* Opening the file */
	FILE *output;

	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}

	/* Writing the output */
	// Header
	fprintf(output, "# step z density\n");

	// Data
	for (int b = 0; b < N_bins; b++)
		fprintf(output, "  %lf %lf\n", z[b], density[b]);

	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}

int extract_sample(int N_configurations, int *N_carbons, Carbon **carbons, Carbon **sample)
{
	printf("Extracting a sample atom...\n");

	/* Allocating the array */
	if ((*sample = malloc(N_configurations * sizeof(Carbon))) == NULL)
	{
		perror("Allocating an array (sample)");
		goto NOMEM;
	}

	/* Extracting the sample atom */
	// Selecting the sample atom
	int serial = carbons[0][0].atom.serial;
	
	// Transfering the data
	for (int c = 0 ; c < N_configurations ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_configurations);
		for (int a = 0 ; a < N_carbons[c]; a++)
		{
			if (carbons[c][a].atom.serial == serial)
			{
				(*sample)[c] = carbons[c][a];
				break;
			}
		}
	}

	/* Success */
	return 0;

	/* Errors */
	NOMEM:
		return ENOMEM;
}

int write_sample(char *file_name, int N_configurations, int *steps, Carbon *sample)
{
	printf("Writing a sample...\n");

	/* Opening the file */
	FILE *output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		goto IO;
	}

	/* Writing to the file */
	fprintf(output, "# step q\n");
	for (int c = 0; c < N_configurations; c++)
		fprintf(output, "  %d %lf\n", steps[c], sample[c].atom.q);
	
	/* Success */
	fclose(output);

	// Exiting normally
	return 0;

	/* Errors */
	IO:
		return EIO;
}

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
	// The plates separation parameters
	const double sep = 3.5;

	// Computing the layers
	compute_layers(N_configurations, box, N_carbons, &carbons, sep);

	/* Computing the charges histogram */
	// The electrodes parameter: indicating the z-coordinates of the lower electrode
	const double limits[2] = {box[0].z_min, box[0].z_min + (box[0].z_max - box[0].z_min) / 2.};
	Electrode electrodes[2] = {LowerElectrode, UpperElectrode};

	// Computing the electrodes
	compute_electrodes(N_configurations, N_carbons, &carbons, limits, electrodes);

	/* Grouping the different carbons */
	// Selecting the lower electrode
	int *N_lower;
	Carbon **lower;
	if ((errno = select_electrode(N_configurations, N_carbons, carbons, LowerElectrode, &N_lower, &lower)) != 0)
		goto EXTRACT_CARBONS;

	// Selecting the inner layer
	int *N_lower_inner;
	Carbon **lower_inner;
	if ((errno = select_layer(N_configurations, N_lower, lower, Inner, &N_lower_inner, &lower_inner)) != 0)
		goto SELECT_LOWER;

	// Selecting the outer layer
	int *N_lower_outer;
	Carbon **lower_outer;
	if ((errno = select_layer(N_configurations, N_lower, lower, Outer, &N_lower_outer, &lower_outer)) != 0)
		goto SELECT_LOWER_INNER;

	// Selecting the upper electrode
	int *N_upper;
	Carbon **upper;
	if ((errno = select_electrode(N_configurations, N_carbons, carbons, UpperElectrode, &N_upper, &upper)) != 0)
		goto SELECT_LOWER_OUTER;

	// Selecting the inner layer
	int *N_upper_inner;
	Carbon **upper_inner;
	if ((errno = select_layer(N_configurations, N_upper, upper, Inner, &N_upper_inner, &upper_inner)) != 0)
		goto SELECT_UPPER;

	// Selecting the outer layer
	int *N_upper_outer;
	Carbon **upper_outer;
	if ((errno = select_layer(N_configurations, N_upper, upper, Outer, &N_upper_outer, &upper_outer)) != 0)
		goto SELECT_UPPER_INNER;

	/* Computing the charges averages */
	// Computing the group
	Group group;

	// The lower electrode average
	if ((errno = average_carbons(N_configurations, N_lower, lower, Q, &group, "The lower electrode")) != 0)
		goto SELECT_UPPER_OUTER;

	if ((errno = write_average("output/graphite/q_lower.dat", N_configurations, steps, group)) != 0)
		goto AVERAGE_CHARGE;

	free(group.average);
	free(group.N);

	// The inner layer of the lower electrode average
	if ((errno = average_carbons(N_configurations, N_lower_inner, lower_inner, Q, &group, "The inner layer of the lower electrode")) != 0)
		goto SELECT_UPPER_OUTER;

	if ((errno = write_average("output/graphite/q_inner-lower.dat", N_configurations, steps, group)) != 0)
		goto AVERAGE_CHARGE;

	free(group.average);
	free(group.N);

	// The outer layer of the lower electrode average
	if ((errno = average_carbons(N_configurations, N_lower_outer, lower_outer, Q, &group, "The outer layer of the lower electrode")) != 0)
		goto SELECT_UPPER_OUTER;

	if ((errno = write_average("output/graphite/q_outer-lower.dat", N_configurations, steps, group)) != 0)
		goto AVERAGE_CHARGE;

	free(group.average);
	free(group.N);

	// The upper electrode average
	if ((errno = average_carbons(N_configurations, N_upper, upper, Q, &group, "The upper electrode")) != 0)
		goto SELECT_UPPER_OUTER;

	if ((errno = write_average("output/graphite/q_upper.dat", N_configurations, steps, group)) != 0)
		goto AVERAGE_CHARGE;

	free(group.average);
	free(group.N);

	// The inner layer of the upper electrode average
	if ((errno = average_carbons(N_configurations, N_upper_inner, upper_inner, Q, &group, "The inner layer of the upper electrode")) != 0)
		goto SELECT_UPPER_OUTER;

	if ((errno = write_average("output/graphite/q_inner-upper.dat", N_configurations, steps, group)) != 0)
		goto AVERAGE_CHARGE;

	free(group.average);
	free(group.N);

	// The outer layer of the upper electrode average
	if ((errno = average_carbons(N_configurations, N_upper_outer, upper_outer, Q, &group, "The outer layer of the upper electrode")) != 0)
		goto SELECT_UPPER_OUTER;

	if ((errno = write_average("output/graphite/q_outer-upper.dat", N_configurations, steps, group)) != 0)
		goto AVERAGE_CHARGE;
	
	/* Computing the average distance between the electrodes' layers' centers of mass */
	Group g1, g2;

	// On the upper electrode
	if ((errno = average_carbons(N_configurations, N_upper_outer, upper_outer, Atom_Z, &g1, "The z-difference between inner and outer layers")) != 0)
		goto AVERAGE_CHARGE;
	
	if ((errno = average_carbons(N_configurations, N_upper_inner, upper_inner, Atom_Z, &g2, "The average z of the inner layer of the upper electrode")) != 0)
		goto AVERAGE_Z1;
	
	group_diff(N_configurations, &g1, g2);

	if ((errno = write_average("output/graphite/dz_upper.dat", N_configurations, steps, g1)) != 0)
		goto AVERAGE_Z2;
	
	free(g1.average), free(g2.average);
	free(g1.N), free(g2.N);

	// On the lower electrode
	if ((errno = average_carbons(N_configurations, N_lower_outer, lower_outer, Atom_Z, &g1, "The z-difference between inner and outer layers")) != 0)
		goto AVERAGE_CHARGE;
	
	if ((errno = average_carbons(N_configurations, N_lower_inner, lower_inner, Atom_Z, &g2, "The average z of the inner layer of the upper electrode")) != 0)
		goto AVERAGE_Z1;
	
	group_diff(N_configurations, &g1, g2);

	if ((errno = write_average("output/graphite/dz_lower.dat", N_configurations, steps, g1)) != 0)
		goto AVERAGE_Z2;

	/* Selecting the sodium */
	int *N_sodium;
	Atom **sodium;
	if ((errno = select_elements(N_configurations, N_atoms, "Na", atoms, &N_sodium, &sodium)) != 0)
		goto AVERAGE_Z2;

	/* Computing the densities */
	double **z, **densities;
	int N_bins = 100;
	if ((errno = compute_density_histograms(N_configurations, N_sodium, sodium, box, &z, &densities, N_bins)) != 0)
		goto SELECT_SODIUM;

	/* Writing the densities */
	if ((errno = write_densities_hist("output/graphite/density_sodium.hist", N_configurations, N_bins, steps, z, densities)) != 0)
		goto DENSITIES;

	double *average_z, *average_density;
	if ((errno = compute_average_density(N_configurations, N_bins, z, densities, &average_z, &average_density)) != 0)
		goto DENSITIES;
	
	if ((errno = write_average_density("output/graphite/average_density_sodium.dat", N_bins, average_z, average_density)) != 0)
		goto AVERAGE_DENSITIES;

	/* Exiting normally */
    exit(EXIT_SUCCESS);

	/* Error handling */
	AVERAGE_DENSITIES:
		free(average_density), free(average_z);
DENSITIES:
	for (int c = 0; c < N_configurations; c++)
		free(densities[c]), free(z[c]);
	free(densities), free(z);
SELECT_SODIUM:
	for (int c = 0; c < N_configurations; c++)
		free(sodium[c]);
	free(sodium), free(N_sodium);
AVERAGE_Z2:
	free(g2.average), free(g2.N);
AVERAGE_Z1:
	free(g1.average), free(g1.N);
AVERAGE_CHARGE:
	free(group.average), free(group.N);
SELECT_UPPER_OUTER:
	for (int c = 0; c < N_configurations; c++)
		free(upper_outer[c]);
	free(upper_outer), free(N_upper_outer);
SELECT_UPPER_INNER:
	for (int c = 0; c < N_configurations; c++)
		free(upper_inner[c]);
	free(upper_inner), free(N_upper_inner);
SELECT_UPPER:
	for (int c = 0; c < N_configurations; c++)
		free(upper[c]);
	free(upper), free(N_upper);
SELECT_LOWER_OUTER:
	for (int c = 0; c < N_configurations; c++)
		free(lower_outer[c]);
	free(lower_outer), free(N_lower_outer);
SELECT_LOWER_INNER:
	for (int c = 0; c < N_configurations; c++)
		free(lower_inner[c]);
	free(lower_inner), free(N_lower_inner);
SELECT_LOWER:
	for (int c = 0; c < N_configurations; c++)
		free(lower[c]);
	free(lower), free(N_lower);
EXTRACT_CARBONS:
	for (int c = 0; c < N_configurations; c++)
		free(carbons[c]);
	free(carbons), free(N_carbons);
READ:
	for (int c = 0; c < N_configurations; c++)
		free(atoms[c]);
	free(atoms);
	free(box), free(N_atoms), free(steps);
ERROR:
	exit(EXIT_FAILURE);
}
