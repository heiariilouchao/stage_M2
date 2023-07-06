#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <argp.h>
#include <math.h>
#define STR_BUFF_LIMIT 256


static char doc[] = "Computing the MSDs of a .lammpstrj configurations file.";

static char args_doc[] = "CONF_FILE";

static struct argp_option options[] =
{
	{
		"start",
		's',
		"START",
		OPTION_ARG_OPTIONAL,
		"The step to start from. Default is 0."
	},
	{0}
};

struct arguments
{
	char *args[1];
	int start;
};

static error_t parse(int key, char *arg, struct argp_state *state)
{
	struct arguments *args = state->input;

	switch (key)
	{
	case 's':
		args->start = atoi(arg);
		if (args->start < 0)
			return EINVAL;
		break;
	case ARGP_KEY_ARG:
		if (state->arg_num >= 1)
			argp_usage(state);
		args->args[state->arg_num] = arg;
		break;
	case ARGP_KEY_END:
		if (state->arg_num < 1)
			argp_usage(state);
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp parser =
{
	options,
	parse,
	args_doc,
	doc
};


error_t read(char *file_name, int timestep, int *N_conf, int **steps, int **N_selection, double ***bounds, int ***indices, double ****pos, double ***charges)
{
	printf("Reading configuration file...\n");


	/* Opening the file */
	FILE *input = fopen(file_name, "r");
	if (input == NULL)
	{
		perror("Opening the configuration file");
		goto IO;
	}


	/* Reading */
	// Initializing
	char str[STR_BUFF_LIMIT];
	const int initial_size = 10;
	int arrays_size = 0, N_atoms;
	*N_conf = 0;

	// Allocating the arrays with the initial size as the number of configurations
	if ((*steps = malloc(initial_size * sizeof(int))) == NULL)
	{
		perror("Allocating an array (steps)");
		goto NOMEM;
	}

	if ((*N_selection = malloc(initial_size * sizeof(int))) == NULL)
	{
		perror("Allocating an array (N_selection)");
		goto NOMEM;
	}

	if ((*bounds = malloc(initial_size * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (bounds)");
		goto NOMEM;
	}

	if ((*indices = malloc(initial_size * sizeof(int *))) == NULL)
	{
		perror("Allocating an array (indices)");
		goto NOMEM;
	}

	if ((*pos = malloc(initial_size * sizeof(double **))) == NULL)
	{
		perror("Allocating an array (indices)");
		goto NOMEM;
	}

	if ((*charges = malloc(initial_size * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (indices)");
		goto NOMEM;
	}

	arrays_size = initial_size;

	// Actually reading
	while (fgets(str, STR_BUFF_LIMIT, input) != NULL)
	{
		if (strcmp(str, "ITEM: TIMESTEP\n") == 0)
		{
			int step;
			if (fscanf(input, "%d\n", &step) != 1)
			{
				perror("Reading the step");
				goto IO;
			}

			while (step < timestep)
			{
				do
				{
					if (fgets(str, STR_BUFF_LIMIT, input) == NULL)
					{
						perror("Skipping the first configurations");
						goto IO;
					}
				}
				while (strcmp(str, "ITEM: TIMESTEP\n") != 0);

				if (fscanf(input, "%d\n", &step) != 1)
				{
					perror("Verifying the timestep");
					goto IO;
				}
			}

			(*N_conf)++;

			if (*N_conf > arrays_size)
				if ((*steps = realloc(*steps, (arrays_size + initial_size) * sizeof(int))) == NULL)
				{
					perror("Resizing an array (steps)");
					goto NOMEM;
				}

			(*steps)[*N_conf - 1] = step;
		}
		else if (strcmp(str, "ITEM: NUMBER OF ATOMS\n") == 0)
		{
			if (fscanf(input, "%d\n", &N_atoms) != 1)
			{
				perror("Reading the number of atoms");
				goto IO;
			}
		}
		else if (strcmp(str, "ITEM: BOX BOUNDS pp pp pp\n") == 0)
		{
			if (*N_conf > arrays_size)
				if ((*bounds = realloc(*bounds, (arrays_size + initial_size) * sizeof(double *))) == NULL)
				{
					perror("Resizing an array (bounds)");
					goto NOMEM;
				}

			if (((*bounds)[*N_conf - 1] = malloc(6 * sizeof(double))) == NULL)
			{
				perror("Allocating an array slot (bounds[])");
				goto NOMEM;
			}

			for (int d = 0 ; d < 3 ; d++)
				if (fscanf(input, "%lf %lf\n", (*bounds)[*N_conf - 1] + 2 * d, (*bounds)[*N_conf - 1] + 2 * d + 1) != 2)
				{
					perror("Reading the box bounds");
					goto IO;
				}
		}
		else if (strcmp(str, "ITEM: ATOMS id element x y z xu yu zu q\n") == 0)
		{
			if (*N_conf > arrays_size)
			{
				if ((*N_selection = realloc(*N_selection, (arrays_size + initial_size) * sizeof(int))) == NULL)
				{
					perror("Resizing an array (N_selection)");
					goto NOMEM;
				}

				if ((*indices = realloc(*indices, (arrays_size + initial_size) * sizeof(int **))) == NULL)
				{
					perror("Resizing an array (indices)");
					goto NOMEM;
				}

				if ((*pos = realloc(*pos, (arrays_size + initial_size) * sizeof(double ***))) == NULL)
				{
					perror("Resizing an array (pos)");
					goto NOMEM;
				}

				if ((*charges = realloc(*charges, (arrays_size + initial_size) * sizeof(double **))) == NULL)
				{
					perror("Resizing an array (charges)");
					goto NOMEM;
				}

				arrays_size += initial_size;
			}

			(*N_selection)[*N_conf - 1] = 0;

			if (((*indices)[*N_conf - 1] = malloc(N_atoms * sizeof(int))) == NULL)
			{
				perror("Allocating an array slot (indices[])");
				goto NOMEM;
			}

			if (((*pos)[*N_conf - 1] = malloc(N_atoms * sizeof(double *))) == NULL)
			{
				perror("Allocating an array slot (pos[])");
				goto NOMEM;
			}

			if (((*charges)[*N_conf - 1] = malloc(N_atoms * sizeof(double))) == NULL)
			{
				perror("Allocating an array slot (charges[])");
				goto NOMEM;
			}

			for (int a = 0 ; a < N_atoms ; a++)
			{
				int index;
				
				if (fscanf(input, "%d %s", &index, str) != 2)
				{
					perror("Reading an atom index and element");
					goto IO;
				}

				if (strncmp(str, "C", 1) != 0)
				{
					if (fgets(str, STR_BUFF_LIMIT, input) == NULL)
					{
						perror("Dumping a line");
						goto IO;
					}

					continue;
				}

				(*N_selection)[*N_conf - 1]++;
				(*indices)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1] = index;

				if (((*pos)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1] = malloc(3 * sizeof(double))) == NULL)
				{
					perror("Allocating an array slot (pos[][])");
					goto NOMEM;
				}

				if (fscanf(input, "%lf %lf %lf %*f %*f %*f %lf",
						   (*pos)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1],
						   (*pos)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1] + 1,
						   (*pos)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1] + 2,
						   (*charges)[*N_conf - 1] + (*N_selection)[*N_conf - 1] - 1) != 4)
				{
					perror("Reading the position and charge");
					goto IO;
				}
			}

			if (((*indices)[*N_conf - 1] = realloc((*indices)[*N_conf - 1], (*N_selection)[*N_conf - 1] * sizeof(int))) == NULL)
			{
				perror("Resizing an array slot (indices[])");
				goto NOMEM;
			}

			if (((*pos)[*N_conf - 1] = realloc((*pos)[*N_conf - 1], (*N_selection)[*N_conf - 1] * sizeof(double *))) == NULL)
			{
				perror("Resizing an array slot (pos[])");
				goto NOMEM;
			}

			if (((*charges)[*N_conf - 1] = realloc((*charges)[*N_conf - 1], (*N_selection)[*N_conf - 1] * sizeof(double))) == NULL)
			{
				perror("Resizing an array slot (charges[])");
			}
		}
		else
		{
			printf("%s\n", str);
			perror("Position in file lost");
			goto IO;
		}
	}


	/* Success */
	// Closing the file
	fclose(input);

	// Exiting normally
	return 0;


	/* Errors */
	IO: return EIO;
	NOMEM: return ENOMEM;
}


error_t write(char *file_name, int N_conf, int *N_selection, int *steps, double **bounds, int **indices, double ***pos, double **charges, int **N_bonds, int ***bonds, bool **layers, bool **electrodes)
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
			fprintf(output, "%d %s", indices[c][a], "C");
			for (int d = 0 ; d < 3 ; d++)
				fprintf(output, " %lf", pos[c][a][d]);
			fprintf(output, " %lf %d", charges[c][a], N_bonds[c][a]);
			for (int b = 0 ; b < N_bonds[c][a] ; b++)
				fprintf(output, " %d", indices[c][bonds[c][a][b]]);
			for (int b = N_bonds[c][a] ; b < 4 ; b++)
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


error_t compute_bonds(int N_conf, int *N_selection, double **bounds, int **indices, double ***pos, double **charges, int ***N_bonds, int ****bonds, const double R, const double delta)
{
	printf("Computing the bonds...\n");


	/* Allocating the arrays */
	if ((*N_bonds = malloc(N_conf * sizeof(int *))) == NULL)
	{
		perror("Allocating an array (N_bonds)");
		return ENOMEM;
	}
	
	if ((*bonds = malloc(N_conf * sizeof(int **))) == NULL)
	{
		perror("Allocating an array (bonds)");
		goto NBONDS;
	}

	for (int c = 0; c < N_conf; c++)
	{
		if (((*N_bonds)[c] = calloc(N_selection[c], sizeof(int))) == NULL)
		{
			perror("Allocating an array slot (N_bonds[])");
			goto BONDS;
		}

		if (((*bonds)[c] = malloc(N_selection[c] * sizeof(int *))) == NULL)
		{
			perror("Allocating an array slot (bonds[])");
			goto BONDS;
		}

		for (int a = 0; a < N_selection[c]; a++)
			if (((*bonds)[c][a] = malloc(4 * sizeof(int))) == NULL)
			{
				perror("Allocating an array slot (bonds[][])");
				goto BONDS;
			}
	}

	/* Actually computing the bonds */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int i = 0 ; i < N_selection[c] ; i++)
		{
			for (int j = i + 1 ; j < N_selection[c] ; j++)
			{
				double r2 = 0.;

				// Applying the PBC
				for (int d = 0 ; d < 3 ; d++)
				{
					double length = bounds[c][2 * d + 1] - bounds[c][2 * d];
					double diff = pos[c][j][d] - pos[c][i][d];

					if (diff < - length / 2.)
						diff += length;
					else if (length / 2. < diff)
						diff -= length;
					
					r2 += diff * diff;
				}

				// if ((R - delta) * (R - delta) <= r2 && r2 <= (R + delta) * (R + delta))
				if (r2 <= (R * R)) // The bond condition
				{
					(*N_bonds)[c][i]++;
					(*bonds)[c][i][(*N_bonds)[c][i] - 1] = j;

					(*N_bonds)[c][j]++;
					(*bonds)[c][j][(*N_bonds)[c][j] - 1] = i;
				}
			}

			// Resizing the array to the actual number of bonds
			if (((*bonds)[c][i] = realloc((*bonds)[c][i], (*N_bonds)[c][i] * sizeof(int))) == NULL)
			{
				perror("Resizing an array slot (bonds[][])");
				goto BONDS;
			}
		}
	}


	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	BONDS: free(bonds);
	NBONDS: free(N_bonds);
	return ENOMEM;
}


error_t compute_layers(int N_conf, int *N_selection, double **bounds, double ***pos, bool ***layers, const double sep, const double delta)
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
				double diff = pos[c][j][2] - pos[c][i][2];

				// Applying the PBC
				double length = bounds[c][2 * 2 + 1] - bounds[c][2 * 2];
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;

				if (!above)
					above = (bool) (sep - delta <= diff) && (diff <= sep + delta);
				if (!under)
					under = (bool) (-sep - delta <= diff) && (diff <= -sep + delta);
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


error_t compute_electrodes(int N_conf, int *N_selection, double **bounds, double ***pos, bool ***electrodes, const double limits[2])
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
			(*electrodes)[c][a] = (limits[0] <= pos[c][a][2] && pos[c][a][2] <= limits[1]);
	}
	
	
	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	ELECTRODES: free(electrodes);
	return ENOMEM;
}


error_t compute_histograms(int N_conf, int *N_selection, bool **layers, bool **electrodes, double **charges, int **N_bonds, double ***histograms, int ***N_types)
{
	printf("Computing the histograms...\n");


	/* Allocating the arrays */
	if ((*histograms = malloc(N_conf * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (histograms)");
		return ENOMEM;
	}

	if ((*N_types = malloc(N_conf * sizeof(int *))) == NULL)
	{
		perror("Allocating an array (N_types)");
		goto HIST;
	}

	for (int c = 0 ; c < N_conf ; c++)
	{
		if (((*histograms)[c] = calloc(5, sizeof(double))) == NULL)	// TODO: number of categories
		{
			perror("Allocating an array slot (histograms[])");
			goto TYPES;
		}

		if (((*N_types)[c] = malloc(5 * sizeof(int))) == NULL)	// TODO: number of categories
		{
			perror("Allocating an array slot (N_types)");
			goto TYPES;
		}
	}


	/* Computing the histograms */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0 ; a < N_selection[c] ; a++)
		{
			int t;

			// TODO: define classification
			if (N_bonds[c][a] == 2)
				t = 4;
			else if (layers[c][a] == 0 && electrodes[c][a] == 0)
				t = 0;
			else if (layers[c][a] == 0 && electrodes[c][a] == 1)
				t = 1;
			else if (layers[c][a] == 1 && electrodes[c][a] == 0)
				t = 2;
			else if (layers[c][a] == 1 && electrodes[c][a] == 1)
				t = 3;
			
			(*histograms)[c][t] += charges[c][a];
			(*N_types)[c][t]++;
		}

		for (int t = 0 ; t < 5 ; t++)	// TODO: number of categories
			(*histograms)[c][t] /= (*N_types)[c][t];
	}


	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	TYPES: free(N_types);
	HIST: free(histograms);
	return ENOMEM;
}


error_t write_hist(char *file_name, int N_conf, int **N_types, int *timestep, double **histograms)
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
	fprintf(output, "# step");
	for (int t = 0 ; t < 5 ; t++)
		fprintf(output, " q%d N%d", t + 1, t + 1);
	fprintf(output, "\n");

	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		fprintf(output, " %d", timestep[c]);
		for (int t = 0 ; t < 5 ; t++)	// TODO: number of categories
			fprintf(output, " %lf %d", histograms[c][t], N_types[c][t]);
		fprintf(output, "\n");
	}


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}


error_t compute_samples(int N_conf, int *N_selection, double **charges, int **N_bonds, bool **layers, bool **electrodes, double ***samples)
{
	/* Allocating the array */
	if ((*samples = malloc(N_conf * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (samples)");
		return ENOMEM;
	}
	
	for (int c = 0 ; c < N_conf ; c++)
	{
		if (((*samples)[c] = malloc(5 * sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (samples[])");
			goto SAMPLES;
		}
	}


	/* Extracting the atoms indices */
	// We assume the groups are static
	int indices[5] = {-1, -1, -1, -1, -1};
	for (int a = 0 ; a < N_selection[0] && (indices[0] == -1 || indices[1] == -1 || indices[2] == -1 || indices[3] == -1 || indices[4] == -1) ; a++)
	{
		if (indices[4] == -1 && N_bonds[0][a] == 2)
			indices[4] = a;
		else if (indices[0] == -1 && (layers[0][a] == 0 && electrodes[0][a] == 0))
			indices[0] = a;
		else if (indices[1] == -1 && (layers[0][a] == 0 && electrodes[0][a] == 1))
			indices[1] = a;
		else if (indices[2] == -1 && (layers[0][a] == 1 && electrodes[0][a] == 0))
			indices[2] = a;
		else if (indices[3] == -1 && (layers[0][a] == 1 && electrodes[0][a] == 1))
			indices[3] = a;
	}


	/* Extracting the charges from the samples */
	for (int c = 0 ; c < N_conf ; c++)
		for (int g = 0 ; g < 5 ; g++)
			(*samples)[c][g] = charges[c][indices[g]];


	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	SAMPLES: free(samples);
	return ENOMEM;
}


error_t write_samples(char *file_name, int N_conf, int *timestep, double **samples)
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
	fprintf(output, "# step");
	for (int t = 0 ; t < 5 ; t++)
		fprintf(output, " q%d", t + 1);
	fprintf(output, "\n");

	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		fprintf(output, " %d", timestep[c]);
		for (int t = 0 ; t < 5 ; t++)	// TODO: number of categories
			fprintf(output, " %lf", samples[c][t]);
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
	/* Error code */
	error_t err;


	/* Parsing the arguments */
	struct arguments arguments;

	// Options' default values
	arguments.start = 0;

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
	{
		perror("Parsing");
		exit(EXIT_FAILURE);
	}


	/* Reading the configurations */
	// Declaring the arrays
	int N_conf, N_atoms, *steps, *N_selection, **indices;
	double ***pos, **bounds, **charges;
	bool *select;

	// Reading the file
	if ((errno = read(arguments.args[0], arguments.start, &N_conf, &steps, &N_selection, &bounds, &indices, &pos, &charges)) != 0)
		goto READ;


	/* Computing the bonds */
	// The bonds parameters
	const double R = 1.7, delta = .5;

	// Declaring the arrays
	int **N_bonds, ***bonds;

	// Computing the bonds
	if ((errno = compute_bonds(N_conf, N_selection, bounds, indices, pos, charges, &N_bonds, &bonds, R, delta)) != 0)
		goto READ;


	/* Computing the layers */
	// The plates separation parameters
	const double sep = 3.2, dz = .5;

	// Declaring the array
	bool **layers;

	// Computing the layers
	if ((errno = compute_layers(N_conf, N_selection, bounds, pos, &layers, sep, dz)) != 0)
		goto BONDS;


	/* Computing the charges histogram */
	// The electrodes parameter
	// const double limit = (bounds[0][5] - bounds[0][4]) / 2.;
	const double limits[2] = {26.0, 35.0};

	// Declaring the array
	bool **electrodes;

	// Computing the electrodes
	if ((errno = compute_electrodes(N_conf, N_selection, bounds, pos, &electrodes, limits)) != 0)
		goto LAYERS;

	
	/* Computing the histograms */
	// Declaring the arrays
	double **histograms;
	int **N_types;

	// Computing the histograms
	if ((errno = compute_histograms(N_conf, N_selection, layers, electrodes, charges, N_bonds, &histograms, &N_types)) != 0)
		goto ELECTRODES;
	
	// Writing the histograms
	if ((errno = write_hist("graphite.hist", N_conf, N_types, steps, histograms)) != 0)
		goto HIST;
	

	/* Computing the samples */
	// Declaring the array
	double **samples;

	// Computing the samples
	if ((errno = compute_samples(N_conf, N_selection, charges, N_bonds, layers, electrodes, &samples)) != 0)
		goto HIST;
	
	// Writing the samples
	if ((errno = write_samples("samples.log", N_conf, steps, samples)) != 0)
		goto SAMPLES;


	/* Writing the output */
	if ((errno = write("graphite.lammpstrj", N_conf, N_selection, steps, bounds, indices, pos, charges, N_bonds, bonds, layers, electrodes)) != 0)
		goto SAMPLES;


	/* Exiting normally */
	exit(EXIT_SUCCESS);


	/* Error handling */
	SAMPLES: free(samples);
	HIST: free(histograms), free(N_types);
	ELECTRODES: free(electrodes);
	LAYERS: free(layers);
	BONDS: free(N_bonds), free(bonds);
	READ: free(select), free(pos);
	exit(EXIT_FAILURE);
}