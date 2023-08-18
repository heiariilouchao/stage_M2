# ifndef CARBONS_H
# define CARBONS_H

# include "../utils/utils.h"

typedef enum Layer
{
	Inner,
	Outer
} Layer;


typedef enum Electrode
{
	LowerElectrode,
	UpperElectrode
} Electrode;


typedef struct Carbon
{
	Atom atom;
	Layer layer;
	Electrode electrode;
} Carbon;

int extract_carbons(int N_configurations, int *N_atoms, Atom **atoms, Box *box, int **N_carbons, Carbon ***carbons);

void compute_layers(int N_conf, Box *box, int *N_carbons, Carbon ***carbons, const double sep);

void compute_electrodes(int N_conf, int *N_carbons, Carbon ***carbons, const double limits[2], Electrode electrode[2]);

int select_electrode(int N_configurations, int *N_carbons, Carbon **carbons, Electrode electrode, int **N_selected, Carbon ***selected);

int select_layer(int N_configurations, int *N_carbons, Carbon **carbons, Layer layer, int **N_selected, Carbon ***selected);

int convert_carbons(int N_configurations, int *N_carbons, Carbon **carbons, Atom ***atoms);

int average_carbons(int N_configurations, int *N_carbons, Carbon **carbons, AtomAttribute attribute, Group *group, char *description);

# endif