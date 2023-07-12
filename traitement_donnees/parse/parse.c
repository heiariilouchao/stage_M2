# include <argp.h>
# include <stdlib.h>
# include <string.h>

# include "parse.h"


struct argp_option options[] =
{
	{
		"start",
		's',
		"START",
		0,
		"The step to start from. Default is 0."
	},
	{
		"bins",
		'n',
		"NBINS",
		0,
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


void set_default_options(Arguments *args)
{
	args->start = 0;
	args->N_bins = 100;
	args->labels = NULL;
	args->elements = NULL;
	args->N_elements = 0;
	args->N_pairs = 0;
}


error_t parse_elements(char *arg, char sep, struct argp_state *state)
{
	printf("Parsing the elements...\n");


	Arguments *args = state->input;
	int n = strlen(arg);

	/* Copying the labels */
	// Checking the labels
	if (arg[0] == sep || arg[n - 1] == sep)
	{
		perror("Invalid labels");
		return EINVAL;
	}
	
	// Actually copying the string
	if ((args->labels = strndup(arg, n)) == NULL)
	{
		perror("Copying a string (args.labels)");
		return ENOMEM;
	}
	

	/* Initializing the other variables */
	args->N_elements = 1;

	if ((args->elements = malloc(sizeof(char *))) == NULL)
	{
		perror("Allocating an array (args.elements)");
		return ENOMEM;
	}
	

	/* Parsing the labels */
	int start = 0;
	for (int c = 0 ; c < n ; c++)
	{
		if (args->labels[c] == sep)
		{
			// Copying an element
			if (((args->elements)[args->N_elements - 1] = strndup(args->labels + start, c - start)) == NULL)
			{
				perror("Copying a string (args.elements[])");
				return ENOMEM;
			}

			(args->N_elements)++;
			start = c + 1;

			// Resizing the array for the next element
			if ((args->elements = realloc(args->elements, args->N_elements * sizeof(char *))) == NULL)
			{
				perror("Resizing an array (args.elements)");
				return ENOMEM;
			}
		}
	}

	// Copying the last element
	if (((args->elements)[args->N_elements - 1] = strndup(args->labels + start, n - start)) == NULL)
	{
		perror("Copying a string (args.elements[])");
		return ENOMEM;
	}

	(args->N_pairs) = (int) args->N_elements * (args->N_elements + 1) / 2;

	printf("Elements informations:\n\tlabels: %s\n\telements: %d", args->labels, args->N_elements);
	for (int e = 0 ; e < args->N_elements ; e++)
		printf(" %s", args->elements[e]);
	printf("\n");


	/* Success */
	return 0;
}


error_t parse(int key, char *arg, struct argp_state *state)
{
	Arguments *args = state->input;
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
		if ((err = parse_elements(arg, ',', state)) != 0)
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


struct argp parser =
{
	options,
	parse,
	args_doc,
	doc
};
