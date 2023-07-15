# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <argp.h>

# include "../parse/parse.h"
# include "../read/read.h"
# include "../bonds/bonds.h"


char doc[] = "Processing and extracting the graphite informations.";

char args_doc[] = "CONF_FILE";


int write(char *file_name, int N_conf, int *N_selection, int *steps, double **bounds, Atom **atoms, bool **layers, bool **electrodes)
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
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);

		fprintf(output, "ITEM: TIMESTEP\n%d\n", steps[c]);

		fprintf(output, "ITEM: NUMBER OF ATOMS\n%d\n", N_selection[c]);

		fprintf(output, "ITEM: BOX BOUNDS pp pp pp\n");
		for (int d = 0 ; d < 3 ; d++)
			fprintf(output, "%lf %lf\n", bounds[c][2 * d], bounds[c][2 * d + 1]);
		
		fprintf(output, "ITEM: ATOMS id element x y z q nb id1 id2 id3 id4 in electrode\n");
		for (int a = 0 ; a < N_selection[c] ; a++)
		{
			fprintf(output, "%d %s", atoms[c][a].serial, "C");
			fprintf(output, " %lf %lf %lf", atoms[c][a].x, atoms[c][a].y, atoms[c][a].z);
			fprintf(output, " %lf %d", atoms[c][a].q, atoms[c][a].N_bonds);
			for (int b = 0 ; b < atoms[c][a].N_bonds ; b++)
				fprintf(output, " %d", atoms[c][atoms[c][a].bonded[b]].serial);
			for (int b = atoms[c][a].N_bonds ; b < 4 ; b++)
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


int compute_layers(int N_conf, int *N_selection, double **bounds, Atom **atoms, bool ***layers, const double sep)
{
	printf("Computing the layers...\n");

	/* Allocating the arrays*/
	if ((*layers = malloc(N_conf * sizeof(bool *))) == NULL)
	{
		perror("Allocating the array (layers)");
		return ENOMEM;
	}

	for (int c = 0; c < N_conf; c++)
		if (((*layers)[c] = calloc(N_selection[c], sizeof(bool))) == NULL)
		{
			perror("Allocating the array slot (layers[c])");
			goto LAYERS;
		}
	

	/* Actually computing the layers */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int i = 0 ; i < N_selection[c] ; i++)
		{
			bool above = false, under = false;

			for (int j = 0; j < N_selection[c] && !(above && under); j++)
			{
				double diff = atoms[c][j].z - atoms[c][i].z;

				// Applying the PBC
				double length = bounds[c][5] - bounds[c][4];
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;

				// TODO: check the conditions !!!
				if (!above)
					above = (bool) (sep / 2. <= diff && diff <= sep);
				if (!under)
					under = (bool) (- sep <= diff && diff <= - sep / 2.);
			}

			(*layers)[c][i] = above && under;
		}
	}
	
	
	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	LAYERS: free(layers);
	return ENOMEM;
}


int compute_electrodes(int N_conf, int *N_selection, double **bounds, Atom **atoms, bool ***electrodes, const double limits[2])
{
	printf("Computing the bonds...\n");


	/* Allocating the electrodes array */
	if ((*electrodes = malloc(N_conf * sizeof(bool *))) == NULL)
	{
		perror("Allocating an array (electrodes)");
		return ENOMEM;
	}

	for (int c = 0 ; c < N_conf ; c++)
		if (((*electrodes)[c] = malloc(N_selection[c] * sizeof(bool))) == NULL)
		{
			perror("Allocating an array slot (electrodes[])");
			goto ELECTRODES;
		}
	


	/* Computing the electrodes */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0 ; a < N_selection[c] ; a++)
			(*electrodes)[c][a] = (limits[0] <= atoms[c][a].z && atoms[c][a].z <= limits[1]);
	}
	
	
	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	ELECTRODES: free(electrodes);
	return ENOMEM;
}


int compute_classification(int N_conf, int *N_selection, bool **layers, bool **electrodes, Atom ***atoms, int *N_groups, Group ***groups, char ***group_descriptions)
{
	printf("Computing the groups...\n");


	/* Initializing the groups */
	// TODO: input the number of groups
	*N_groups = 5;

	// Initializing the array
	if ((*groups = malloc(N_conf * sizeof(Group *))) == NULL)
	{
		perror("Allocating an array (groups)");
		return ENOMEM;
	}

	for (int c = 0 ; c < N_conf ; c++)
	{
		if (((*groups)[c] = malloc(*N_groups * sizeof(Group))) == NULL)
		{
			perror("Allocating an array slot (groups[])");
			goto GROUPS;
		}

		for (int g = 0 ; g < *N_groups ; g++)
		{
			(*groups)[c][g].N = 0;
			(*groups)[c][g].average = 0.;
		}
	}


	/* Initializing the group descriptions */
	// TODO: input group descriptions
	const char * descriptions[] =
	{
		"Outer layer of negative electrode",
		"Outer layer of positive electrode",
		"Inner layer of negative electrode",
		"Inner layer of positive electrode",
		"Bonded to two atoms"
	};

	// Copying the descriptions
	if ((*group_descriptions = malloc(*N_groups * sizeof(char *))) == NULL)
	{
		perror("Allocating an array (group_descriptions)");
		goto GROUPS;
	}

	for (int g = 0 ; g < *N_groups ; g++)
	{
		if (((*group_descriptions)[g] = malloc(STR_GROUP_DESC_LIMIT * sizeof(char))) == NULL)
		{
			perror("Allocating an array slot (group_descriptions[])");
			goto DESCRIPTIONS;
		}

		if (strcmp(strncpy((*group_descriptions)[g], descriptions[g], STR_GROUP_DESC_LIMIT - 1), descriptions[g]) != 0)
		{
			printf("strings:\n%s\n%s\n", descriptions[g], (*group_descriptions)[g]);
			perror("Copying a string (group_descriptions[])");
			goto DESCRIPTIONS;
		}
	}


	/* Computing the classification */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0 ; a < N_selection[c] ; a++)
		{
			int t;

			// TODO: define the classification
			if ((*atoms)[c][a].N_bonds == 2)
				t = 4;
			else
			{
				if (layers[c][a] == 0 && electrodes[c][a] == 0)
					t = 0;
				else if (layers[c][a] == 0 && electrodes[c][a] == 1)
					t = 1;
				else if (layers[c][a] == 1 && electrodes[c][a] == 0)
					t = 2;
				else if (layers[c][a] == 1 && electrodes[c][a] == 1)
					t = 3;
			}
			
			(*atoms)[c][a].group = t;
			(*groups)[c][t].N++;
		}
	}


	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	DESCRIPTIONS: free(group_descriptions);
	GROUPS: free(groups);
	return ENOMEM;
}


int compute_charges_histograms(int N_conf, int *N_selection, Atom **atoms, bool **layers, bool **electrodes, int N_groups, Group ***groups)
{
	printf("Computing the histograms...\n");

	
	/* Computing the histograms */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0 ; a < N_selection[c] ; a++)
			(*groups)[c][atoms[c][a].group].average += atoms[c][a].q;

		for (int t = 0 ; t < N_groups ; t++)
			(*groups)[c][t].average /= (*groups)[c][t].N;
	}


	/* Success */
	// Exiting normally
	return 0;
}


int compute_density_histograms(int N_conf, int *N_selection, Atom **atoms, double **bounds, double ***z, double ***densities, const int N_bins)
{
	printf("Computing the density histogram...\n");

	double delta = (bounds[0][5] - bounds[0][4]) / N_bins;


	/* Allocating the array */
	// int *N;
	// if ((N = calloc(N_conf, sizeof(int))) == NULL)
	// {
	// 	perror("Allocating an array (N)");
	// 	return ENOMEM;
	// }

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

	for (int c = 0 ; c < N_conf ; c++)
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

		double z0 = bounds[c][4];
		for (int b = 0 ; b < N_bins ; b++)
			(*z)[c][b] = z0 + delta * b + delta / 2.;
	}

	


	/* Computing the densities */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0 ; a < N_selection[c] ; a++)
		{
			int bin = (atoms[c][a].z - bounds[c][4]) / delta;
			if (bin < N_bins)
			{
				// N[c]++;
				(*densities)[c][bin]++;
			}
		}

		for (int b = 0 ; b < N_bins ; b++)
			(*densities)[c][b] /= (bounds[c][1] - bounds[c][0]) * (bounds[c][3] - bounds[c][2]) * delta;
	}


	/* Success */
	return 0;


	/* Errors */
	DENSITIES: free(densities);
	Z: free(z);
	// N: free(N);
	return ENOMEM;
}


int write_charges_hist(char *file_name, int N_conf, int *timestep, int N_groups, char **group_descriptions, Group **groups)
{
	printf("Writing the histograms...\n");


	/* Opening the file */
	FILE* output;

	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing the output */
	// Groups informations
	fprintf(output, "# Groups informations:\n");
	for (int g = 0 ; g < N_groups ; g++)
		fprintf(output, "#   %d: %s\n", g + 1, group_descriptions[g]);
	
	// Header
	fprintf(output, "# step");
	for (int t = 0 ; t < 5 ; t++)
		fprintf(output, " q%d N%d", t + 1, t + 1);
	fprintf(output, "\n");

	// Data
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		fprintf(output, " %d", timestep[c]);
		for (int t = 0 ; t < N_groups ; t++)
			fprintf(output, " %lf %d", groups[c][t].average, groups[c][t].N);
		fprintf(output, "\n");
	}


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}


int write_densities_hist(char *file_name, int N_conf, int N_bins, int *timestep, double **z, double **densities)
{
	printf("Writing the densities histograms...\n");


	/* Opening the file */
	FILE* output;

	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing the output */
	// Header
	fprintf(output, "# step z density\n");

	// Data
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int b = 0 ; b < N_bins ; b++)
			fprintf(output, "  %d %lf %lf\n", timestep[c], z[c][b], densities[c][b]);
		fprintf(output, "\n\n");
	}


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}


int compute_samples(int N_conf, int *N_selection, Atom **atoms, int N_groups, Atom ***samples)
{
	/* Allocating the array */
	int *indices;
	if ((indices = malloc(N_groups * sizeof(int))) == NULL)
	{
		perror("Allocating an array (indices)");
		return ENOMEM;
	}

	for (int t = 0 ; t < N_groups ; t++)
		indices[t] = -1;

	if ((*samples = malloc(N_conf * sizeof(Atom *))) == NULL)
	{
		perror("Allocating an array (samples)");
		goto INDICES;
	}
	
	for (int c = 0 ; c < N_conf ; c++)
	{
		if (((*samples)[c] = malloc(N_groups * sizeof(Atom))) == NULL)
		{
			perror("Allocating an array slot (samples[])");
			goto SAMPLES;
		}
	}


	/* Extracting the atoms indices */
	// We assume the groups are static
	bool done = false;
	for (int a = 0 ; a < N_selection[0]; a++)
	{
		done = true;
		for (int t = 0 ; t < N_groups ; t++)
			done = done && (indices[t] != -1);
		
		if (done)
			break;

		int t = atoms[0][a].group;
		if (indices[t] == -1)
			indices[t] = a;
	}


	/* Extracting the charges from the samples */
	for (int c = 0 ; c < N_conf ; c++)
		for (int t = 0 ; t < N_groups ; t++)
			(*samples)[c][t] = atoms[c][indices[t]];


	/* Success */
	free(indices);

	// Exiting normally
	return 0;


	/* Errors */
	INDICES: free(indices);
	SAMPLES: free(samples);
	return ENOMEM;
}


error_t write_samples(char *file_name, int N_conf, int *timestep, int N_groups, Atom **samples)
{
	printf("Writing the samples...\n");


	/* Opening the file */
	FILE* output;

	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing the output */
	fprintf(output, "# Samples serial indices:\n");
	for (int t = 0 ; t < N_groups ; t++)
		fprintf(output, "#   %d: %d\n", t + 1, samples[0][t].serial);

	fprintf(output, "# step");
	for (int t = 0 ; t < N_groups ; t++)
		fprintf(output, " q%d", t + 1);
	fprintf(output, "\n");

	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		fprintf(output, " %d", timestep[c]);
		for (int t = 0 ; t < N_groups ; t++)
			fprintf(output, " %lf", samples[c][t].q);
		fprintf(output, "\n");
	}


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}


int main(int argc, char **argv)
{
	/* Parsing the arguments */
	Arguments arguments;

	// Options' default values
	set_default_options(&arguments);
	arguments.labels = "C";
	arguments.elements = malloc(sizeof(char *));
	arguments.elements[0] = "C";
	arguments.N_elements = 1;


	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
	{
		perror("Parsing");
		exit(EXIT_FAILURE);
	}


	/* Reading the carbons */
	// Declaring the arrays
	int N_conf_carbons, *steps_carbons, *N_carbons;
	Atom **carbons;
	double **bounds_carbons;

	// Reading the file
	if ((errno = read_trajectory(&arguments, &N_conf_carbons, &steps_carbons, &N_carbons, &bounds_carbons, &carbons)) != 0)
		goto ERROR;


	/* Computing the bonds */
	// The bonds parameters
	const double R = 1.7;

	// Computing the bonds
	if ((errno = compute_cutoff_bonds(N_conf_carbons, N_carbons, bounds_carbons, &carbons, R)) != 0)
		goto READ_CARBONS;


	/* Computing the layers */
	// The plates separation parameters
	const double sep = 3.5;

	// Declaring the array
	bool **layers;

	// Computing the layers
	if ((errno = compute_layers(N_conf_carbons, N_carbons, bounds_carbons, carbons, &layers, sep)) != 0)
		goto READ_CARBONS;


	/* Computing the charges histogram */
	// The electrodes parameter
	const double limits[2] = {26.0, 35.0};

	// Declaring the array
	bool **electrodes;

	// Computing the electrodes
	if ((errno = compute_electrodes(N_conf_carbons, N_carbons, bounds_carbons, carbons, &electrodes, limits)) != 0)
		goto LAYERS;

	
	/* Computing the histograms */
	// Declaring the arrays
	int N_groups;
	Group **groups;
	char **group_descriptions;

	// Computing the histograms
	if ((errno = compute_classification(N_conf_carbons, N_carbons, layers, electrodes, &carbons, &N_groups, &groups, &group_descriptions)) != 0)
		goto ELECTRODES;

	if ((errno = compute_charges_histograms(N_conf_carbons, N_carbons, carbons, layers, electrodes, N_groups, &groups)) != 0)
		goto GROUPS;
	
	// Writing the histograms
	if ((errno = write_charges_hist("output/graphite.hist", N_conf_carbons, steps_carbons, N_groups, group_descriptions, groups)) != 0)
		goto GROUPS;
	

	/* Computing the samples */
	// Declaring the array
	Atom **samples;

	// Computing the samples
	if ((errno = compute_samples(N_conf_carbons, N_carbons, carbons, N_groups, &samples)) != 0)
		goto GROUPS;
	
	// Writing the samples
	if ((errno = write_samples("output/samples.log", N_conf_carbons, steps_carbons, N_groups, samples)) != 0)
		goto SAMPLES;


	/* Writing the output */
	if ((errno = write("output/graphite.lammpstrj", N_conf_carbons, N_carbons, steps_carbons, bounds_carbons, carbons, layers, electrodes)) != 0)
		goto SAMPLES;

	
	/* Reading the sodium */
	int N_conf_sodium, *steps_sodium, *N_sodium;
	Atom **sodium;
	double **bounds_sodium;

	arguments.labels = "Na";
	arguments.elements[0] = "Na";
	if ((errno = read_trajectory(&arguments, &N_conf_sodium, &steps_sodium, &N_sodium, &bounds_sodium, &sodium)) != 0)
		goto SAMPLES;
	

	/* Computing the densities */
	double **z, **densities;
	int N_bins = 50;
	if ((errno = compute_density_histograms(N_conf_sodium, N_sodium, sodium, bounds_sodium, &z, &densities, N_bins)) != 0)
		goto READ_SODIUM;
	

	/* Writing the densities */
	if ((errno = write_densities_hist("output/sodium.hist", N_conf_sodium, N_bins, steps_sodium, z, densities)) != 0)
		goto DENSITIES;


	/* Exiting normally */
	exit(EXIT_SUCCESS);


	/* Error handling */
	DENSITIES: free(densities), free(z);
	READ_SODIUM: free(bounds_sodium), free(sodium), free(N_sodium), free(steps_sodium);
	SAMPLES: free(samples);
	GROUPS: free(group_descriptions), free(groups);
	ELECTRODES: free(electrodes);
	LAYERS: free(layers);
	READ_CARBONS: free(bounds_carbons), free(carbons), free(N_carbons), free(steps_carbons);
	ERROR: exit(EXIT_FAILURE);
}