# include <argp.h>
# include <stdlib.h>
# include <string.h>

# include "parse.h"
# include "../utils/utils.h"


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
	args->labels = NULL;
	args->elements = NULL;
	args->N_elements = 0;
}


int parse_elements(char *labels, int *N_elements, char ***elements)
{
    printf("Parsing the elements...\n");
    
    
    /* Allocating the array */
	if ((*elements = malloc(sizeof(char **))) == NULL)
    {
        perror("Allocating an array (elements)");
        return ENOMEM;
    }

	char *element = strtok(labels, ",");
	*N_elements = 0;
	while (element)
	{
		if (((*elements)[*N_elements] = malloc(STR_ELEMENT_LIMIT * sizeof(char))) == NULL)
		{
			perror("Allocating an array slot (elements[])");
			goto NOMEM;
		}

		memcpy((*elements)[*N_elements], element, STR_ELEMENT_LIMIT - 1);

		(*N_elements)++;
		element = strtok(NULL, ",");
	}

    /* Success */
    // Prompting the informations
    printf("Elements informations:\n\tN_elements: %d\n\telements:", *N_elements);
    for (int e = 0 ; e < *N_elements ; e++)
        printf(" %s", (*elements)[e]);
    printf("\n");
    
    // Exiting normally
    return 0;


    /* Errors */
    NOMEM: free(elements);
    return ENOMEM;
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
	case 'l':
		if ((err = parse_elements(arg, &((*args).N_elements), &((*args).elements))) != 0)
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
