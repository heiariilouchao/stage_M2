# ifndef PARSER_H
# define PARSER_H

# include <argp.h>

extern char doc[];

extern char args_doc[];

extern struct argp_option options[];

struct arguments
{
	char *args[1];
	int start, N_bins;
	char *labels, **elements;
	int N_elements, N_pairs;
};


extern error_t parse_elements(char *arg, char sep, struct arguments **args);

extern error_t parse(int key, char *arg, struct argp_state *state);


extern struct argp parser;


# endif
