#include <stdio.h>
#include <stdlib.h>
#include <errno.h>


#include "bonds.h"
#include "read.h"


int compute_cutoff_bonds(int N_conf, int *N_selection, double **bounds, Atom ***atoms, const double R)
{
	printf("Computing the bonds with cutoff...\n");
    printf("\tR = %lf\n", R);


	/* Initializing the values */
	for (int c = 0 ; c < N_conf ; c++)
		for (int a = 0 ; a < N_selection[c] ; a++)
		{
			(*atoms)[c][a].N_bonds = 0;
			if (((*atoms)[c][a].bonded = malloc(4 * sizeof(int))) == NULL)
			{
				perror("Allocating an array (atoms[][].bonded)");
				return ENOMEM;
			}
		}

	/* Actually computing the bonds */
	for (int c = 0 ; c < N_conf ; c++)
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int i = 0 ; i < N_selection[c] ; i++)
		{
			for (int j = i + 1 ; j < N_selection[c] ; j++)
			{
				double r2 = 0.;
				double length, diff;

				length = bounds[c][1] - bounds[c][0];
				diff = (*atoms)[c][j].x - (*atoms)[c][i].x;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				length = bounds[c][3] - bounds[c][2];
				diff = (*atoms)[c][j].y - (*atoms)[c][i].y;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				length = bounds[c][5] - bounds[c][4];
				diff = (*atoms)[c][j].z - (*atoms)[c][i].z;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				if (r2 <= (R * R)) // The bond condition
				{
					(*atoms)[c][i].N_bonds++;
					(*atoms)[c][i].bonded[(*atoms)[c][i].N_bonds - 1] = j;

					(*atoms)[c][j].N_bonds++;
					(*atoms)[c][j].bonded[(*atoms)[c][j].N_bonds - 1] = i;
				}
			}

			// Resizing the array to the actual number of bonds
			// if (((*bonds)[c][i] = realloc((*bonds)[c][i], (*N_bonds)[c][i] * sizeof(int))) == NULL)
			if (((*atoms)[c][i].bonded = realloc((*atoms)[c][i].bonded, (*atoms)[c][i].N_bonds * sizeof(int))) == NULL)
			{
				perror("Resizing an array slot (atoms[][].bonded)");
				goto BONDS;
			}
		}
	}


	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	BONDS:
	for (int c = 0 ; c < N_conf ; c++)
		for (int a = 0 ; a < N_selection[c] ; a++)
			free((*atoms)[c][a].bonded);
	return ENOMEM;
}
