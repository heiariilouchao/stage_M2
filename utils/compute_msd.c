# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <argp.h>
# include <math.h>
# define STR_BUFF_LIMIT 256
# define STR_LABELS_LIMIT 3
# define STR_N_LABELS_LIMIT 4


static char doc[] = "Computing the MSDs of a .lammpstrj configurations file.";

static char args_doc[] = "CONF_FILE";

static struct argp_option options[] =
{
	{
		"steps",
		's',
		"MAX",
		OPTION_ARG_OPTIONAL,
		"The maximum number of configurations to process. Default is 0, meaning all."
	},
	{
		"start",
		't',
		"START",
		OPTION_ARG_OPTIONAL,
		"The step to start from. Default is 0."
	},
	{
		"labels",
		'l',
		"LABEL1,LABEL2,...",
		OPTION_ARG_OPTIONAL,
		"The labels of the atoms to select for the computation. Default is all."
	},
	{0}
};

struct arguments
{
	char *args[1], *labels;
	int max, start, N_labels;
};

static error_t parse(int key, char *arg, struct argp_state *state)
{
	struct arguments *args = state->input;
	
	switch (key)
	{
		case 's':
			args->max = atoi(arg);
			if (args->max < 0)
				return EINVAL;
			break;
		case 'k':
			args->start = atoi(arg);
			if (args->start < 0)
				return EINVAL;
			break;
		case 'l':
			args->N_labels = 1;
			args->labels = (char *) realloc(args->labels, (strlen(arg) + 1) * sizeof(char));
			if (strcpy(args->labels, arg) == NULL)
				return EINVAL;

			for (int c = 0 ; c < strlen(args->labels) ; c++)
			{
				if (args->labels[c] == ',')
					(args->N_labels)++;
			}
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


error_t read(char *file_name, int timestep, int *N_conf, int *N_atoms, double ****upos, int N_labels, char *labels, bool **select)
{
	/* Opening the file */
	FILE* conf_file_ptr = fopen(file_name, "r");
	
	// Error handling
	if (conf_file_ptr == NULL)
	{
		perror("Opening the configuration file to read the parameters");
		goto IO;
	}

	
	char dump[STR_BUFF_LIMIT], str_value[STR_BUFF_LIMIT];
	int error, step = 0;
	long file_pos;

	/* Skipping the first configurations */
	while (step < timestep) // Based on the current timestep
	{
		while (strcmp(str_value, "ITEM: TIMESTEP\n") != 0)  // Based on the header
		{
			if (fgets(str_value, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
			{
				perror("Skipping (flag)");
				goto IO;
			}
		}
		if (fscanf(conf_file_ptr, "%d\n", &step) != 1)
		{
			perror("Skipping (step)");
			goto IO;
		}
	}


	/* Pre-reading the file */
	// Searching for the number of atoms
	while (strcmp(str_value, "ITEM: NUMBER OF ATOMS\n") != 0) // Based on the header
		if (fgets(str_value, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
		{
			perror("Pre-reading number of atoms (flag)");
			goto IO;
		}
	if (fscanf(conf_file_ptr, "%d\n", N_atoms) != 1)
	{
		perror("Pre-reading (N_atoms)");
		goto IO;
	}


	// Searching for the start of the configuration in the file
	while (strncmp(str_value, "ITEM: ATOMS", 11) != 0)  // Based on the header
	{
		if (fgets(str_value, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
		{
			perror("Pre-reading start of file (flag)");
			goto IO;
		}
	}
	if ((file_pos = ftell(conf_file_ptr)) == -1L)
	{
		perror("Pre-reading start of file (file_pos)");
		goto IO;
	}

	// Selecting the atoms
	*select = (bool *) malloc(*N_atoms * sizeof(bool));
	if (labels != NULL)
	{
		for (int a = 0 ; a < *N_atoms ; a++)
		{
			int index;
			if (fscanf(conf_file_ptr, "%d %s", &index, str_value) != 2)	// Reading first two fields
			{
				perror("Pre-reading atom selection (index)");
				goto IO;
			}
			if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)	// Dumping rest of line
			{
				perror("Pre-reading atom selection (dump)");
				goto IO;
			}
			
			index--;
			(*select)[index] = false;
			if (strstr(labels, str_value) != NULL)
				(*select)[index] = true;
		}

		if (fseek(conf_file_ptr, file_pos, SEEK_SET) != 0)	// Repositioning in the file
		{
			perror("Resetting the positions to start of file (select)");
			goto IO;
		}
	}
	else
		for (int a = 0 ; a < *N_atoms ; a++)
			(*select)[a] = true;


	// Counting the number of configurations to process
	*N_conf = 1;
	while (fgets(str_value, STR_BUFF_LIMIT, conf_file_ptr) != NULL)
		if (strncmp(str_value, "ITEM: ATOMS", 11) == 0)	// Based on file position
			(*N_conf)++;
	

	/* Allocating the array */
	if ((*upos = (double ***) malloc(*N_conf * sizeof(double **))) == NULL)
	{
		perror("Allocating the array (upos)");
		goto NOMEM;
	}


	/* Returning to the start of the file */
	if (fseek(conf_file_ptr, file_pos, SEEK_SET) != 0)
	{
		perror("Resetting the positions to start of file");
		goto IO;
	}


	/* Actually reading the configurations */
	// Going through all configurations
	for (int c = 0 ; c < *N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, *N_conf);
		// Allocating the configuration slot in the array
		if (((*upos)[c] = (double **) malloc(*N_atoms * sizeof(double *))) == NULL)
		{
			perror("Reading the configurations (upos[c])");
			goto NOMEM;
		}

		// Going through all atoms
		for (int a = 0 ; a < *N_atoms ; a++)
		{
			// Reading the index of the current atom
			int index;
			if (fscanf(conf_file_ptr, "%d", &index) != 1)
			{
				perror("Reading the configurations (index)");
				goto IO;
			}
			index--;
				
			// Allocating the atom slot in the array
			if (((*upos)[c][index] = (double *) malloc(3 * sizeof(double))) == NULL)
			{
				perror("Reading the configurations (upos[c][index])");
				goto NOMEM;
			}
			
			// Reading the position
			if (fscanf(conf_file_ptr, 
				"%*s %*f %*f %*f %lf %lf %lf %*f\n",
				(*upos)[c][index], (*upos)[c][index] + 1, (*upos)[c][index] + 2) != 3)
			{
				perror("Reading the configurations (upos)");
				goto IO;
			}
		}

		// Going to the next configuration
		if (c != *N_conf - 1)
		{
			do
			{
				if (fgets(str_value, STR_BUFF_LIMIT, conf_file_ptr) == NULL)	// Dumping a line
				{
					perror("Reading the configurations (flag)");
					goto IO;
				}
			}
			while (strncmp(str_value, "ITEM: ATOMS", 11) != 0); // Based on the header
		}
	}

	// Closing the file
	fclose(conf_file_ptr);
	return 0;

	// Error handling
	IO: return EIO;
	NOMEM: return ENOMEM;
}


error_t write(char *file_name, double *MSD, int N_conf)
{
	FILE* output = fopen(file_name, "w");

	if (output == NULL)
	{
		perror("Opening the output file");
		return EIO;
	}

	fprintf(output, "# tau\tMSD\n");
	fprintf(output, " 0\t0.0\n");

	for (int tau = 1 ; tau < N_conf ; tau++)
	{
		printf("tau: %d / %d\r", tau, N_conf);
		fprintf(output, " %d\t%lf\n", tau, MSD[tau]);
	}
	
	fclose(output);
	return 0;
}


int main(int argc, char **argv)
{
	/* Error code */
	error_t err;


	/* Parsing the arguments */
	struct arguments arguments;

	// Options' default values
	arguments.max = 0;
	arguments.start = 0;
	arguments.labels = NULL;
	arguments.N_labels = 0;

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
	{
		perror("Parsing");
		exit(EXIT_FAILURE);
	}


	/* Reading the parameters */
	printf("Reading configuration file...\n");
	int N_conf, N_atoms;
	double ***upos;
	bool *select;
	if ((errno = read(arguments.args[0], arguments.start, &N_conf, &N_atoms, 
			  &upos, arguments.N_labels, arguments.labels, &select)) != 0)
		goto READ;

	
	/* Computing the MSD */
	// Allocating the array
	printf("Computing MSD...\n");
	double *MSD;
	if ((MSD = (double *) calloc(N_conf - 1, sizeof(double))) == NULL)
	{
		perror("Allocating the MSD array");
		goto READ;
	}

	for (int tau = 1 ; tau < N_conf ; tau++)
	{
		printf("tau: %d / %d\r", tau, N_conf);
		
		for (int c = 0 ; c < N_conf - tau ; c++)
			for (int a = 0 ; a < N_atoms ; a++)
				if (select[a])
				{
					double dist = 0.;
					for (int d = 0 ; d < 2 ; d++)
						dist += (upos[c + tau][a][d] - upos[c][a][d]) * (upos[c + tau][a][d] - upos[c][a][d]);
					MSD[tau] += sqrt(dist);
				}

		MSD[tau] = 1. / (N_conf - tau) * 1. / N_atoms * MSD[tau];
	}


	/* Writing the full MSD */
	printf("Writing output...\n");
	if ((errno = write("msd.dat", MSD, N_conf)) != 0)
	{
		perror("Writing the MSD");
		goto WRITE;
	}


	/* Exiting normally */
	exit(EXIT_SUCCESS);


	/* Error handling */
	WRITE: free(MSD);
	READ: free(select), free(upos);
	exit(EXIT_FAILURE);
}