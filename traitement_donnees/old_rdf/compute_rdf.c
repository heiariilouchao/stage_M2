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
	{
		"bins",
		'n',
		"NBINS",
		OPTION_ARG_OPTIONAL,
		"The number of bins to compute the RDF. Default is 100."
	},
	{
		"labels",
		'l',
		"LABEL1,LABEL2,...",
		0,
		"The labels of the atoms to select for the computation."
	},
	{0}
};

struct arguments
{
	char *args[1];
	int start, N_bins;
	char *labels, **elements;
	int N_elements, N_pairs;
};


error_t parse_elements(char *arg, char sep, struct arguments **args)
{
	int n = strlen(arg);

	/* Copying the labels */
	// Checking the labels
	if (arg[0] == sep || arg[n - 1] == sep)
	{
		perror("Invalid labels");
		return EINVAL;
	}
	
	// Actually copying the string
	if (((*args)->labels = strndup(arg, n)) == NULL)
	{
		perror("Copying a string (args.labels)");
		return ENOMEM;
	}
	

	/* Initializing the other variables */
	(*args)->N_elements = 1;

	if (((*args)->elements = malloc(sizeof(char *))) == NULL)
	{
		perror("Allocating an array (args.elements)");
		return ENOMEM;
	}
	

	/* Parsing the labels */
	int start = 0;
	for (int c = 0 ; c < n ; c++)
	{
		if ((*args)->labels[c] == sep)
		{
			// Copying an element
			if ((((*args)->elements)[(*args)->N_elements - 1] = strndup((*args)->labels + start, c - start)) == NULL)
			{
				perror("Copying a string (args.elements[])");
				return ENOMEM;
			}

			((*args)->N_elements)++;
			start = c + 1;

			// Resizing the array for the next element
			if (((*args)->elements = realloc((*args)->elements, (*args)->N_elements * sizeof(char *))) == NULL)
			{
				perror("Resizing an array (args.elements)");
				return ENOMEM;
			}
		}
	}

	// Copying the last element
	if ((((*args)->elements)[(*args)->N_elements - 1] = strndup((*args)->labels + start, n - start)) == NULL)
	{
		perror("Copying a string (args.elements[])");
		return ENOMEM;
	}

	((*args)->N_pairs) = (int) (*args)->N_elements * ((*args)->N_elements + 1) / 2;


	/* Success */
	return 0;
}


static error_t parse(int key, char *arg, struct argp_state *state)
{
	struct arguments *args = state->input;
	error_t err = 0;

	switch (key)
	{
	case 's':
		if (arg < 0)
			return EINVAL;
		args->start = atoi(arg);
		break;
	case 'n':
		if (arg < 0)
			return EINVAL;
		args->N_bins = atoi(arg);
		break;
	case 'l':
		if ((err = parse_elements(arg, ',', &args)) != 0)
			return err;
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


// error_t read(char *file_name, int timestep, int N_elements, char *labels, char **elements, int *N_conf, int ***N_selection, double ***bounds, double *****pos)
// {
// 	printf("Reading configuration file...\n");


// 	/* Opening the file */
// 	FILE *input = fopen(file_name, "r");
// 	if (input == NULL)
// 	{
// 		perror("Opening the configuration file");
// 		goto IO;
// 	}


// 	/* Reading */
// 	// Initializing
// 	char str[STR_BUFF_LIMIT], hdr[STR_BUFF_LIMIT];
// 	const int initial_size = 10;
// 	int arrays_size = 0, N_atoms;
// 	*N_conf = 0;

// 	// Allocating the arrays with the initial size as the number of configurations
// 	if ((*N_selection = malloc(initial_size * sizeof(int *))) == NULL)
// 	{
// 		perror("Allocating an array (N_selection)");
// 		goto NOMEM;
// 	}

// 	if ((*bounds = malloc(initial_size * sizeof(double *))) == NULL)
// 	{
// 		perror("Allocating an array (bounds)");
// 		goto NOMEM;
// 	}

// 	if ((*pos = malloc(initial_size * sizeof(double ***))) == NULL)
// 	{
// 		perror("Allocating an array (indices)");
// 		goto NOMEM;
// 	}

// 	arrays_size = initial_size;

// 	// Actually reading
// 	while (fgets(hdr, STR_BUFF_LIMIT, input) != NULL)
// 	{
// 		if (strcmp(hdr, "ITEM: TIMESTEP\n") == 0)
// 		{
// 			int step;
// 			if (fscanf(input, "%d\n", &step) != 1)
// 			{
// 				perror("Reading the step");
// 				goto IO;
// 			}

// 			while (step < timestep)
// 			{
// 				do
// 				{
// 					if (fgets(hdr, STR_BUFF_LIMIT, input) == NULL)
// 					{
// 						perror("Skipping the first configurations");
// 						goto IO;
// 					}
// 				}
// 				while (strcmp(hdr, "ITEM: TIMESTEP\n") != 0);

// 				if (fscanf(input, "%d\n", &step) != 1)
// 				{
// 					perror("Verifying the timestep");
// 					goto IO;
// 				}
// 			}

// 			(*N_conf)++;
// 		}
// 		else if (strcmp(hdr, "ITEM: NUMBER OF ATOMS\n") == 0)
// 		{
// 			if (fscanf(input, "%d\n", &N_atoms) != 1)
// 			{
// 				perror("Reading the number of atoms");
// 				goto IO;
// 			}
// 		}
// 		else if (strcmp(hdr, "ITEM: BOX BOUNDS pp pp pp\n") == 0)
// 		{
// 			if (*N_conf > arrays_size)
// 				if ((*bounds = realloc(*bounds, (arrays_size + initial_size) * sizeof(double *))) == NULL)
// 				{
// 					perror("Resizing an array (bounds)");
// 					goto NOMEM;
// 				}

// 			if (((*bounds)[*N_conf - 1] = malloc(6 * sizeof(double))) == NULL)
// 			{
// 				perror("Allocating an array slot (bounds[])");
// 				goto NOMEM;
// 			}

// 			for (int d = 0 ; d < 3 ; d++)
// 				if (fscanf(input, "%lf %lf\n", (*bounds)[*N_conf - 1] + 2 * d, (*bounds)[*N_conf - 1] + 2 * d + 1) != 2)
// 				{
// 					perror("Reading the box bounds");
// 					goto IO;
// 				}
// 		}
// 		else if (strncmp(hdr, "ITEM: ATOMS id element", 22) == 0)
// 		{
// 			if (*N_conf > arrays_size)
// 			{
// 				if ((*N_selection = realloc(*N_selection, (arrays_size + initial_size) * sizeof(int *))) == NULL)
// 				{
// 					perror("Resizing an array (N_selection)");
// 					goto NOMEM;
// 				}

// 				if ((*pos = realloc(*pos, (arrays_size + initial_size) * sizeof(double ***))) == NULL)
// 				{
// 					perror("Resizing an array (pos)");
// 					goto NOMEM;
// 				}

// 				arrays_size += initial_size;
// 			}

// 			if (((*N_selection)[*N_conf - 1] = calloc(N_elements, sizeof(int))) == NULL)
// 			{
// 				perror("Allocating an array slot (N_selection[][])");
// 				goto NOMEM;
// 			}

// 			if (((*pos)[*N_conf - 1] = malloc(N_elements * sizeof(double **))) == NULL)
// 			{
// 				perror("Allocating an array slot (pos[])");
// 				goto NOMEM;
// 			}

// 			for (int e = 0 ; e < N_elements ; e++)
// 				if (((*pos)[*N_conf - 1][e] = malloc(N_atoms * sizeof(double *))) == NULL)
// 				{
// 					perror("Allocating an array slot (pos[][])");
// 					goto NOMEM;
// 				}

// 			for (int a = 0 ; a < N_atoms ; a++)
// 			{
// 				int type;
				
// 				if (fscanf(input, "%*d %s", str) != 1)
// 				{
// 					perror("Reading an atom element");
// 					goto IO;
// 				}

// 				if (strstr(labels, str) == NULL)
// 				{
// 					if (fgets(str, STR_BUFF_LIMIT, input) == NULL)
// 					{
// 						perror("Dumping a line");
// 						goto IO;
// 					}

// 					continue;
// 				}

// 				for (int e = 0 ; e < N_elements ; e++)
// 					if (strcmp(str, elements[e]) == 0)
// 					{
// 						type = e;
// 						break;
// 					}

// 				(*N_selection)[*N_conf - 1][type]++;

// 				if (((*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1] = malloc(3 * sizeof(double))) == NULL)
// 				{
// 					perror("Allocating an array slot (pos[][])");
// 					goto NOMEM;
// 				}

// 				if (strcmp(hdr, "ITEM: ATOMS id element x y z xu yu zu\n") == 0)
// 					if (fscanf(input, "%lf %lf %lf %*f %*f %*f\n",
// 							(*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1],
// 							(*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1] + 1,
// 							(*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1] + 2) != 3)
// 					{
// 						perror("Reading the position");
// 						goto IO;
// 					}
// 				else if (strcmp(hdr, "ITEM: ATOMS id element x y z xu yu zu q\n") == 0)
// 					if (fscanf(input, "%lf %lf %lf %*f %*f %*f %*f\n",
// 							(*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1],
// 							(*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1] + 1,
// 							(*pos)[*N_conf - 1][type][(*N_selection)[*N_conf - 1][type] - 1] + 2) != 3)
// 					{
// 						perror("Reading the position");
// 						goto IO;
// 					}
// 				else
// 				{
// 					perror("Unknown dumping format");
// 					goto IO;
// 				}
// 			}

// 			for (int e = 0 ; e < N_elements ; e++)
// 				if (((*pos)[*N_conf - 1][e] = realloc((*pos)[*N_conf - 1][e], (*N_selection)[*N_conf - 1][e] * sizeof(double *))) == NULL)
// 				{
// 					perror("Resizing an array slot (pos[][])");
// 					goto NOMEM;
// 				}
// 		}
// 		else
// 		{
// 			printf("%s\n", hdr);
// 			perror("Position in file lost");
// 			goto IO;
// 		}
// 	}


// 	/* Success */
// 	// Closing the file
// 	fclose(input);

// 	// Exiting normally
// 	return 0;


// 	/* Errors */
// 	IO: return EIO;
// 	NOMEM: return ENOMEM;
// }


// int transpose(int N_elements, int e1, int e2)
// {
// 	int sum = 0;
// 	for (int e = 0 ; e < N_elements ; e++)
// 		sum += N_elements - e;
	
// 	return sum + e2 - e1;
// }


error_t compute_rdf(int N_conf, int N_elements, int **N_selection, double **bounds, double ****pos, int N_bins, double cutoff, int N_pairs, double **r, double ***RDF)
{
	double delta = cutoff / N_bins;

	/* Allocating the arrays */
	int **hist;

	if (((*r) = malloc(N_bins * sizeof(double))) == NULL)
	{
		perror("Allocating an array (r)");
		return ENOMEM;
	}

	if (((hist = malloc(N_pairs * sizeof(double *)))) == NULL)
	{
		perror("Allocating an array (hist)");
		goto R;
	}

	if (((*RDF) = malloc(N_pairs * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (RDF)");
		goto HIST;
	}

	for (int p = 0 ; p < N_pairs ; p++)
	{
		if ((hist[p] = calloc(N_bins, sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (hist[])");
			goto RDF;
		}

		if (((*RDF)[p] = calloc(N_bins, sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (RDF[])");
			goto RDF;
		}
	}

	// Initializing the distance range
	for (int b = 0 ; b < N_bins ; b++)
		(*r)[b] = b * delta;


	/* Incrementing the histograms */
	printf("Incrementing the histograms...\n");

	for (int c = 0 ; c < N_conf ; c++)
		for (int e1 = 0, p = 0 ; e1 < N_elements ; e1++)
			for (int e2 = e1 ; e2 < N_elements ; e2++, p++)
				for (int a1 = 0 ; a1 < N_selection[c][e1] ; a1++)
					for (int a2 = (e1 == e2) ? a1 + 1 : 0 ; a2 < N_selection[c][e2] ; a2++)
					{
						double r2 = 0.;
						for (int d = 0 ; d < 3 ; d++)
						{
							double length = bounds[c][2 * d + 1] - bounds[c][2 * d];
							double diff = pos[c][e2][a2][d] - pos[c][e1][a1][d];
							if (diff < - length / 2.)
								diff += length;
							else if (length / 2. < diff)
								diff -= length;
							
							r2 += diff * diff;
						}

						int bin = (int) (sqrt(r2) / delta);
						if (bin < N_bins)
							hist[p][bin] += (e1 == e2) ? 2 : 1;
					}
	

	/* Computing the RDFs */
	printf("Normalizing the RDFs...\n");

	double V = 1.;
	for (int d = 0 ; d < 3 ; d++)
		V *= bounds[0][2 *d + 1] - bounds[0][2 * d];
	double coeff = 1. / (4. / 3. * M_PI / V);

	for (int e1 = 0, p = 0 ; e1 < N_elements ; e1++)
		for (int e2 = e1 ; e2 < N_elements ; e2++, p++)
		{
			int N = (e1 == e2) ? N_selection[0][e1] * (N_selection[0][e1] - 1) : N_selection[0][e1] * N_selection[0][e2];
			for (int b = 0 ; b < N_bins ; b++)
				(*RDF)[p][b] = (double) coeff * hist[p][b] / N_conf / (N * (pow((*r)[b] + delta, 3) - pow((*r)[b], 3)));
		}
	

	/* Success */
	return 0;


	/* Errors */
	RDF: free(RDF);
	HIST: free(hist);
	R: free(r);
	return ENOMEM;
}


error_t write(char *file_name, int N_elements, char **elements, int N_bins, int N_pairs, double *r, double **RDF)
{
	printf("Writing the output...\n");

	/* Opening the file */
	FILE* output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing */
	fprintf(output, "# r");
	for (int e1 = 0 ; e1 < N_elements ; e1++)
		for (int e2 = e1 ; e2 < N_elements ; e2++)
			fprintf(output, " g_%s%s", elements[e1], elements[e2]);
	fprintf(output, "\n");


	for (int b = 0 ; b < N_bins ; b++)
	{
		fprintf(output, "  %lf", r[b]);
		for (int p = 0 ; p < N_pairs ; p++)
			fprintf(output, " %lf", RDF[p][b]);
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
	arguments.N_bins = 100;

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
		exit(EXIT_FAILURE);


	/* Reading the configurations */
	// Declaring the arrays
	int N_conf, **N_selection;
	double ****pos, **bounds;

	// Reading the file
	if ((errno = read(arguments.args[0], arguments.start, arguments.N_elements, arguments.labels, arguments.elements, &N_conf, &N_selection, &bounds, &pos)) != 0)
		exit(EXIT_FAILURE);
	
	

	/* Computing the RDFs */
	double *r, **RDF;
	double cutoff = 10.0;
	if ((errno = compute_rdf(N_conf, arguments.N_elements, N_selection, bounds, pos, arguments.N_bins, cutoff, arguments.N_pairs, &r, &RDF)) != 0)
		goto READ;
	

	/* Writing the output */
	if ((errno = write("rdf.dat", arguments.N_elements, arguments.elements, arguments.N_bins, arguments.N_pairs, r, RDF)) != 0)
		goto RDF;


	/* Exiting normally */
	exit(EXIT_SUCCESS);


	/* Error handling */
	RDF: free(RDF), free(r);
	READ: free(pos), free(bounds), free(N_selection);
	exit(EXIT_FAILURE);
}
