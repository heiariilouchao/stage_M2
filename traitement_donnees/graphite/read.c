# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>

# include "utils.h"
# include "read.h"


int read_trajectory(char *file_name, int timestep, int *N_conf, int **steps, int **N_selection, double ***bounds, Atom ***atoms)
{
	printf("Reading configuration file...\n");


	/* Opening the file */
	FILE *input;
	if ((input = fopen(file_name, "r")) == NULL)
	{
		perror("Opening the configuration file");
		goto IO;
	}


	/* Reading */
	// Initializing
	char str[STR_BUFF_LIMIT];
	int arrays_size = 0, N_atoms;
	*N_conf = 0;

	// Allocating the arrays with the initial size as the number of configurations
	if ((*steps = malloc(INITIAL_CONF * sizeof(int))) == NULL)
	{
		perror("Allocating an array (steps)");
		goto NOMEM;
	}

	if ((*N_selection = malloc(INITIAL_CONF * sizeof(int))) == NULL)
	{
		perror("Allocating an array (N_selection)");
		goto NOMEM;
	}

	if ((*bounds = malloc(INITIAL_CONF * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (bounds)");
		goto NOMEM;
	}

    if ((*atoms = malloc(INITIAL_CONF * sizeof(Atom *))) == NULL)
	{
		perror("Allocating an array (atoms)");
		goto NOMEM;
	}

	arrays_size = INITIAL_CONF;

	// Actually reading
	while (fgets(str, STR_BUFF_LIMIT, input) != NULL)
	{
		if (strcmp(str, "ITEM: TIMESTEP\n") == 0)
		{
			int step;
			if (fscanf(input, "%d\n", &step) != 1)
			{
				perror("Reading the step");
				goto IO;
			}

			while (step < timestep)
			{
				do
				{
					if (fgets(str, STR_BUFF_LIMIT, input) == NULL)
					{
						perror("Skipping the first configurations");
						goto IO;
					}
				}
				while (strcmp(str, "ITEM: TIMESTEP\n") != 0);

				if (fscanf(input, "%d\n", &step) != 1)
				{
					perror("Verifying the timestep");
					goto IO;
				}
			}

			(*N_conf)++;

			if (*N_conf > arrays_size)
				if ((*steps = realloc(*steps, (arrays_size + INCREMENT_CONF) * sizeof(int))) == NULL)
				{
					perror("Resizing an array (steps)");
					goto NOMEM;
				}

			(*steps)[*N_conf - 1] = step;
		}
		else if (strcmp(str, "ITEM: NUMBER OF ATOMS\n") == 0)
		{
			if (fscanf(input, "%d\n", &N_atoms) != 1)
			{
				perror("Reading the number of atoms");
				goto IO;
			}
		}
		else if (strcmp(str, "ITEM: BOX BOUNDS pp pp pp\n") == 0)
		{
			if (*N_conf > arrays_size)
				if ((*bounds = realloc(*bounds, (arrays_size + INCREMENT_CONF) * sizeof(double *))) == NULL)
				{
					perror("Resizing an array (bounds)");
					goto NOMEM;
				}

			if (((*bounds)[*N_conf - 1] = malloc(6 * sizeof(double))) == NULL)
			{
				perror("Allocating an array slot (bounds[])");
				goto NOMEM;
			}

			for (int d = 0 ; d < 3 ; d++)
				if (fscanf(input, "%lf %lf\n", (*bounds)[*N_conf - 1] + 2 * d, (*bounds)[*N_conf - 1] + 2 * d + 1) != 2)
				{
					perror("Reading the box bounds");
					goto IO;
				}
		}
		else if (strcmp(str, "ITEM: ATOMS id element x y z xu yu zu q\n") == 0)
		{
			if (*N_conf > arrays_size)
			{
				if ((*N_selection = realloc(*N_selection, (arrays_size + INCREMENT_CONF) * sizeof(int))) == NULL)
				{
					perror("Resizing an array (N_selection)");
					goto NOMEM;
				}

                if ((*atoms = realloc(*atoms, (arrays_size + INCREMENT_CONF) * sizeof(Atom *))) == NULL)
                {
					perror("Resizing an array (atoms)");
					goto NOMEM;
				}

				arrays_size += INCREMENT_CONF;
			}

			(*N_selection)[*N_conf - 1] = 0;

            if (((*atoms)[*N_conf - 1] = malloc(N_atoms * sizeof(Atom))) == NULL)
            {
				perror("Allocating an array slot (atoms[])");
				goto NOMEM;
			}

			for (int a = 0 ; a < N_atoms ; a++)
			{
				int index;
				
				if (fscanf(input, "%d %s", &index, str) != 2)
				{
					perror("Reading an atom index and element");
					goto IO;
				}

				int error;
				if ((error = select_atom(input, "C", str)) == -1)
					goto IO;
				else if (error == 1)
					continue;

				(*N_selection)[*N_conf - 1]++;
                (*atoms)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1].serial = index;

				if (fscanf(input, "%lf %lf %lf %*f %*f %*f %lf",
                    &((*atoms)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1].x),
                    &((*atoms)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1].y),
                    &((*atoms)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1].z),
                    &((*atoms)[*N_conf - 1][(*N_selection)[*N_conf - 1] - 1].q)) != 4)
				{
					perror("Reading the position and charge");
					goto IO;
				}
			}

            if (((*atoms)[*N_conf - 1] = realloc((*atoms)[*N_conf - 1], (*N_selection)[*N_conf - 1] * sizeof(Atom))) == NULL)
            {
                perror("Resizing an array slot (atoms[])");
				goto NOMEM;
            }
		}
		else
		{
			printf("%s\n", str);
			perror("Position in file lost");
			goto IO;
		}
	}


	/* Success */
	// Closing the file
	fclose(input);

	// Exiting normally
	return 0;


	/* Errors */
	IO: return EIO;
	NOMEM: return ENOMEM;
}


int select_atom(FILE *file, char *labels, char *element)
{
	if (strstr(labels, element) == NULL)
	{
		if (fgets(element, STR_BUFF_LIMIT, file) == NULL)
		{
			perror("Dumping a line");
			return -1;
		}

		return 1;
	}

	return 0;
}
