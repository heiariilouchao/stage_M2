#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <errno.h>

#include "../utils/utils.h"
#include "../bonds/bonds.h"
#include "carbons.h"

int extract_carbons(int N_configurations, int *N_atoms, Atom **atoms, Box *box, int **N_carbons, Carbon ***carbons)
{
	printf("Extracting the carbons...\n");

	/* Selecting the carbons */
	int *N_carbons_local;
	Atom **carbons_local;

	if (select_elements(N_configurations, N_atoms, "C", atoms, &N_carbons_local, &carbons_local) != 0)
		goto ERROR;

	if (compute_cutoff_bonds(N_configurations, N_carbons_local, box, &carbons_local, 1.7) != 0)
		goto SELECT_ELEMENTS;

	/* Transfering the data */
	if ((*N_carbons = malloc(N_configurations * sizeof(int))) == NULL)
	{
		perror("Allocating an array (N_carbons)");
		goto SELECT_ELEMENTS;
	}

	if ((*carbons = malloc(N_configurations * sizeof(Carbon *))) == NULL)
	{
		perror("Allocating an array (carbons)");
		goto N_CARBONS;
	}

	for (int c = 0; c < N_configurations; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_configurations);

		if (((*carbons)[c] = malloc(N_carbons_local[c] * sizeof(Carbon))) == NULL)
		{
			perror("Allocating an array slot (carbons[])");
			goto CARBONS;
		}

		for (int a = 0; a < N_carbons_local[c]; a++)
			(*carbons)[c][a].atom = carbons_local[c][a];
		(*N_carbons)[c] = N_carbons_local[c];
	}

	/* Success */
	// Freeing the temporary variables
	for (int c = 0; c < N_configurations; c++)
		free(carbons_local[c]);
	free(carbons_local), free(N_carbons_local);

	// Exiting normally
	return 0;

/* Errors */
CARBONS:
	free(carbons);
N_CARBONS:
	free(N_carbons);
SELECT_ELEMENTS:
	for (int c = 0; c < N_configurations; c++)
		free(carbons_local[c]);
	free(carbons_local), free(N_carbons_local);
ERROR:
	return 1;
}

void compute_layers(int N_conf, Box *box, int *N_carbons, Carbon ***carbons, const double sep)
{
	printf("Computing the layers...\n");

	/* Actually computing the layers */
	for (int c = 0; c < N_conf; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int i = 0; i < N_carbons[c]; i++)
		{
			bool above = false, under = false;

			for (int j = 0; j < N_carbons[c] && !(above && under); j++)
			{
				double diff = (*carbons)[c][j].atom.z - (*carbons)[c][i].atom.z;

				// Applying the PBC
				double length = box[c].z_max - box[c].z_min;
				if (diff < -length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;

				if (!above)
					above = (bool)(sep / 2. <= diff && diff <= sep);
				if (!under)
					under = (bool)(-sep <= diff && diff <= -sep / 2.);
			}

			if (above && under)
				(*carbons)[c][i].layer = Inner;
			else
				(*carbons)[c][i].layer = Outer;
		}
	}
}

void compute_electrodes(int N_conf, int *N_carbons, Carbon ***carbons, const double limits[2], Electrode electrode[2])
{
	printf("Computing the bonds...\n");

	/* Computing the electrodes */
	for (int c = 0; c < N_conf; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int a = 0; a < N_carbons[c]; a++)
			if (limits[0] <= (*carbons)[c][a].atom.z && (*carbons)[c][a].atom.z <= limits[1])
				(*carbons)[c][a].electrode = electrode[0];
			else
				(*carbons)[c][a].electrode = electrode[1];
	}
}

int select_electrode(int N_configurations, int *N_carbons, Carbon **carbons, Electrode electrode, int **N_selected, Carbon ***selected)
{
	printf("Selecting an electrode...\n");

	/* Allocating the arrays */
	if ((*N_selected = calloc(N_configurations, sizeof(int))) == NULL)
	{
		perror("Allocating an array (N_selected)");
		goto NOMEM;
	}

	if ((*selected = malloc(N_configurations * sizeof(Carbon *))) == NULL)
	{
		perror("Allocating an array (selected)");
		goto N_SELECTED;
	}

	/* Selecting the carbons */
	for (int c = 0; c < N_configurations; c++)
	{
		(*selected)[c] = NULL;
		for (int a = 0; a < N_carbons[c]; a++)
			if (carbons[c][a].electrode == electrode)
			{
				(*N_selected)[c]++;

				if (((*selected)[c] = realloc((*selected)[c], (*N_selected)[c] * sizeof(Carbon))) == NULL)
				{
					perror("Allocating an array slot (selected[])");
					goto SELECTED;
				}

				(*selected)[c][(*N_selected)[c] - 1] = carbons[c][a];
			}
	}

	/* Success */
	return 0;

/* Errors */
SELECTED:
	free(*selected);
N_SELECTED:
	free(*N_selected);
NOMEM:
	return ENOMEM;
}

int select_layer(int N_configurations, int *N_carbons, Carbon **carbons, Layer layer, int **N_selected, Carbon ***selected)
{
	printf("Selecting a layer...\n");

	/* Allocating the arrays */
	if ((*N_selected = calloc(N_configurations, sizeof(int))) == NULL)
	{
		perror("Allocaitng an array (N_selected)");
		goto NOMEM;
	}

	if ((*selected = malloc(N_configurations * sizeof(Carbon *))) == NULL)
	{
		perror("Allocating an array (selected)");
		goto N_SELECTED;
	}

	/* Selecting the carbons */
	for (int c = 0; c < N_configurations; c++)
	{
		(*selected)[c] = NULL;
		for (int a = 0; a < N_carbons[c]; a++)
			if (carbons[c][a].layer == layer)
			{
				(*N_selected)[c]++;

				if (((*selected)[c] = realloc((*selected)[c], (*N_selected)[c] * sizeof(Carbon))) == NULL)
				{
					perror("Reallocating an array slot (selected[])");
					goto SELECTED;
				}

				(*selected)[c][(*N_selected)[c] - 1] = carbons[c][a];
			}
	}

	/* Prompting the informations */
	char *layer_str;
	switch (layer)
	{
		case Outer:
			layer_str = "outer layer";
			break;
		case Inner:
			layer_str = "inner layer";
			break;
	}
	printf("Selection informations:\n\tSelection: %s\n\tN_selected: %d\n", layer_str, (*N_selected)[0]);

	/* Success */
	return 0;

/* Errors */
SELECTED:
	free(*selected);
N_SELECTED:
	free(*N_selected);
NOMEM:
	return ENOMEM;
}

int convert_carbons(int N_configurations, int *N_carbons, Carbon **carbons, Atom ***atoms)
{
	/* Allocating the array */
	if ((*atoms = malloc(N_configurations * sizeof(Atom *))) == NULL)
	{
		perror("Allocating an array (atoms)");
		goto NOMEM;
	}

	for (int c = 0 ; c < N_configurations ; c++)
	{
		if (((*atoms)[c] = malloc(N_carbons[c] * sizeof(Atom))) == NULL)
		{
			perror("Allocating an array slot (atoms[])");
			goto ATOMS;
		}
		
		for (int a = 0 ; a < N_carbons[c] ; a++)
			(*atoms)[c][a] = carbons[c][a].atom;
	}

	/* Success */
	// Exiting normally
	return 0;

	/* Errors */
	ATOMS:
		free(atoms);
	NOMEM:
		return ENOMEM;
}
