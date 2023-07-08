# include <argp.h>
# include <stdlib.h>

# include "parse.h"


char doc[] = "Extracting, processing informations about supercapacitors electrodes from a file.";

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
	{0}
};


error_t parse(int key, char *arg, struct argp_state *state)
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


struct argp parser =
{
	options,
	parse,
	args_doc,
	doc
};
