# ifndef PARSER_H
# define PARSER_H

# include <argp.h>

extern char doc[];

extern char args_doc[];

extern struct argp_option options[];

typedef struct Arguments
{
	char *args[1];
	int start;
	char *labels;
	int N_elements;
	char **elements;
} Arguments;


extern void set_default_options(Arguments *args);

int parse_elements(char *labels, int *N_elements, char ***elements);

extern error_t parse(int key, char *arg, struct argp_state *state);


extern struct argp parser;

# endif
