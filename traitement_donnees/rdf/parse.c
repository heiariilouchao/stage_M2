# include <argp.h>
# include <stdlib.h>
# include <string.h>

# include "parse.h"


char doc[] = "Computing the MSDs of a .lammpstrj configurations file.";

char args_doc[] = "CONF_FILE";

struct argp_option options[] =
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


error_t parse(int key, char *arg, struct argp_state *state)
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


struct argp parser =
{
	options,
	parse,
	args_doc,
	doc
};
