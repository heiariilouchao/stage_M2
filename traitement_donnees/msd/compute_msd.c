# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <argp.h>
# include <math.h>
# include <omp.h>
# define STR_BUFF_LIMIT 256


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
		"labels",
		'l',
		"LABEL1,LABEL2,...",
		OPTION_ARG_OPTIONAL,
		"The labels of the atoms to select for the computation. Default is all (NULL)."
	},
	{0}
};

struct arguments
{
	char *args[1], *labels;
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
		case 'l':
			args->labels = (char *) realloc(args->labels, (strlen(arg) + 1) * sizeof(char));
			if (strcpy(args->labels, arg) == NULL)
				return EINVAL;
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 1)
				argp_usage(state);
			args->args[state->arg_num] = arg;
			break;
		case ARGP_KEY_END:
			if (state->arg_num <1)
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


error_t read(char *file_name, int timestep, int *N_conf, int **N_selection, char *selection, int ***indices, double ****upos)
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
	char hdr[STR_BUFF_LIMIT], str[STR_BUFF_LIMIT];
	const int initial_size = 10;
	int arrays_size = 0, N_atoms;
	*N_conf = 0;

	// Allocating the arrays with the initial size as the number of configurations
	if ((*N_selection = malloc(initial_size * sizeof(int))) == NULL)
	{
		perror("Allocating an array (N_selection)");
		goto NOMEM;
	}

	if ((*indices = malloc(initial_size * sizeof(int *))) == NULL)
	{
		perror("Allocating an array (indices)");
		goto NOMEM;
	}

	if ((*upos = malloc(initial_size * sizeof(double **))) == NULL)
	{
		perror("Allocating an array (upos)");
		goto NOMEM;
	}

	arrays_size = initial_size;

	// Actually reading
	while (fgets(hdr, STR_BUFF_LIMIT, input) != NULL)
	{
		if (strcmp(hdr, "ITEM: TIMESTEP\n") == 0)
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
					if (fgets(hdr, STR_BUFF_LIMIT, input) == NULL)
					{
						perror("Skipping the first configurations");
						goto IO;
					}
				}
				while (strcmp(hdr, "ITEM: TIMESTEP\n") != 0);

				if (fscanf(input, "%d\n", &step) != 1)
				{
					perror("Verifying the timestep");
					goto IO;
				}
			}

			(*N_conf)++;
		}
		else if (strcmp(hdr, "ITEM: NUMBER OF ATOMS\n") == 0)
		{
			if (fscanf(input, "%d\n", &N_atoms) != 1)
			{
				perror("Reading the number of atoms");
				goto IO;
			}
		}
		else if (strcmp(hdr, "ITEM: BOX BOUNDS pp pp pp\n") == 0)
		{
			for (int d = 0 ; d < 3 ; d++)
				if (fgets(str, STR_BUFF_LIMIT, input) == NULL)
				{
					perror("Dumping the box bounds");
					goto IO;
				}
		}
		// TODO: verify the dumping format
		else if (strncmp(hdr, "ITEM: ATOMS id element", 22) == 0)
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

				if ((*upos = realloc(*upos, (arrays_size + initial_size) * sizeof(double ***))) == NULL)
				{
					perror("Resizing an array (upos)");
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

			if (((*upos)[*N_conf - 1] = malloc(N_atoms * sizeof(double *))) == NULL)
			{
				perror("Allocating an array slot (upos[])");
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

				if (selection != NULL && strstr(selection, str) == NULL)
				{
					if (fgets(str, STR_BUFF_LIMIT, input) == NULL)
					{
						perror("Dumping a line");
						goto IO;
					}

					continue;
				}

				index--;
				(*N_selection)[*N_conf - 1]++;
				(*indices)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1] = index;

				if (((*upos)[*N_conf - 1][index] = malloc(3 * sizeof(double))) == NULL)
				{
					perror("Allocating an array slot (upos[][])");
					goto NOMEM;
				}

				if (strcmp(hdr, "ITEM: ATOMS id element x y z xu yu zu\n") == 0)
					if (fscanf(input, "%*f %*f %*f %lf %lf %lf\n",
							(*upos)[*N_conf - 1][index],
							(*upos)[*N_conf - 1][index] + 1,
							(*upos)[*N_conf - 1][index] + 2) != 3)
					{
						perror("Reading the position and charge");
						goto IO;
					}
				else if (strcmp(hdr, "ITEM: ATOMS id element x y z xu yu zu q\n") == 0)
					if (fscanf(input, "%*f %*f %*f %lf %lf %lf %*f\n",
							(*upos)[*N_conf - 1][index],
							(*upos)[*N_conf - 1][index] + 1,
							(*upos)[*N_conf - 1][index] + 2) != 3)
					{
						perror("Reading the position and charge");
						goto IO;
					}
				else
				{
					perror("Unknown dumping format");
					goto IO;
				}
			}

			if (((*indices)[*N_conf - 1] = realloc((*indices)[*N_conf - 1], (*N_selection)[*N_conf - 1] * sizeof(int))) == NULL)
			{
				perror("Resizing an array slot (indices[])");
				goto NOMEM;
			}

			// upos n'est pas réalloué car on l'utilise comme sparse array en se basant sur les indices des atomes
		}
		else
		{
			printf("%s\n", hdr);
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


error_t compute_msd(int N_conf, int *N_selection, int **indices, double ***upos, double **MSD)
{
	printf("Computing MSD...\n");


	/* Allocating the array */
	if ((*MSD = calloc(N_conf, sizeof(double))) == NULL)
	{
		perror("Allocating an array (MSD)");
		return ENOMEM;
	}


	/* Actually computing the MSD */
	#pragma omp parallel shared(N_conf, N_selection, indices, upos, MSD)
	{
		#pragma omp for schedule(dynamic)
		for (int tau = 1 ; tau < N_conf ; tau++)
		{
			for (int c = 0 ; c < N_conf - tau ; c++)
			{
				double sum = 0.;

				for (int a = 0 ; a < N_selection[c] ; a++)
					for (int d = 0 ; d < 3 ; d++)
						sum += (upos[c + tau][indices[c][a]][d] - upos[c][indices[c][a]][d]) * (upos[c + tau][indices[c][a]][d] - upos[c][indices[c][a]][d]);

				(*MSD)[tau] += (double) sum / N_selection[c];
			}

			(*MSD)[tau] /= (double) (N_conf - tau);
		}
	}
}


error_t write(char *file_name, double *MSD, int N_conf)
{
	printf("Writing the MSD...\n");


	/* Opening the output file */
	FILE* output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening the output file");
		return EIO;
	}


	/* Actually writing */
	// Writing the header
	fprintf(output, "# tau\tMSD\n");

	for (int tau = 0 ; tau < N_conf ; tau++)
		fprintf(output, " %d\t%lf\n", tau, MSD[tau]);
	
	
	/* Success */
	// Closing the output file
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
	arguments.labels = NULL;

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
	{
		perror("Parsing");
		exit(EXIT_FAILURE);
	}


	/* Reading the parameters */
	int N_conf, *N_selection, **indices;
	double ***upos;
	if ((errno = read(arguments.args[0], arguments.start, &N_conf, &N_selection, arguments.labels, &indices, &upos)) != 0)
		goto READ;

	
	/* Computing the MSD */
	// Allocating the array
	double *MSD;
	if ((errno = compute_msd(N_conf, N_selection, indices, upos, &MSD)) != 0)
		goto READ;


	/* Writing the full MSD */
	if ((errno = write("msd.dat", MSD, N_conf)) != 0)
		goto COMPUTE;


	/* Exiting normally */
	exit(EXIT_SUCCESS);


	/* Errors */
	COMPUTE: free(MSD);
	READ: free(indices), free(N_selection), free(upos);
	exit(EXIT_FAILURE);
}
