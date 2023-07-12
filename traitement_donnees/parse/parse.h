# ifndef PARSER_H
# define PARSER_H

# include <argp.h>

extern char doc[];

extern char args_doc[];

extern struct argp_option options[];

typedef struct Arguments
{
	char *args[1];
	int start, N_bins;
	char *labels, **elements;
	int N_elements, N_pairs;
} Arguments;


extern void set_default_options(Arguments *args);

extern error_t parse_elements(char *arg, char sep, struct argp_state *state);

extern error_t parse(int key, char *arg, struct argp_state *state);


extern struct argp parser;


# endif
