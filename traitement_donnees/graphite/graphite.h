# ifndef GRAPHITE_H
# define GRAPHITE_H


# include "../utils/utils.h"


typedef enum Layer
{
	Inner,
	Outer
} Layer;


typedef enum Electrode
{
	Negative,
	Positive
} Electrode;


typedef struct Carbon
{
	Atom atom;
	Layer layer;
	Electrode electrode;
} Carbon;

# endif