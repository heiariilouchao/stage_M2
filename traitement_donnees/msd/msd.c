# include <stdio.h>
# include <stdlib.h>
# include <errno.h>
# include <omp.h>

# include "../utils.h"
# include "../parse/parse.h"
# include "../read/read.h"


char doc[] = "Computing the MSDs of a .lammpstrj configurations file.";

char args_doc[] = "CONF_FILE";


int compute_msd(int N_conf, int *N_selection, Atom **atoms, double **MSD)
{
	printf("Computing the MSD...\n");


	/* Allocating the array */
	if ((*MSD = calloc(N_conf, sizeof(double))) == NULL)
	{
		perror("Allocating an array (MSD)");
		return ENOMEM;
	}


	/* Actually computing the MSD */
	# pragma omp parallel shared(N_conf, N_selection, atoms, MSD)
	{
		# ifdef _OPENMP
		# pragma omp master
		printf("Using %d threads for computation...\n", omp_get_num_threads());
		# endif
		# pragma omp for schedule(dynamic)
		for (int tau = 1 ; tau < N_conf ; tau++)
		{
			printf("tau: %d / %d\r", tau, N_conf - 1);
			for (int c = 0 ; c < N_conf - tau ; c++)
			{
				double sum = 0.;

				for (int i = 0 ; i < N_selection[c + tau] ; i++)
				{
					Atom atom = atoms[c + tau][i];
					for (int j = 0 ; j < N_selection[c] ; j++)
						if (atom.serial == atoms[c][j].serial)
						{
							sum += (atom.xu - atoms[c][j].xu) * (atom.xu - atoms[c][j].xu) + (atom.yu - atoms[c][j].yu) * (atom.yu - atoms[c][j].yu) + (atom.zu - atoms[c][j].zu) * (atom.zu - atoms[c][j].zu);
							break;
						}
				}

				(*MSD)[tau] += (double) sum / N_selection[c];
			}

			(*MSD)[tau] /= (double) (N_conf - tau);
		}
	}


	/* Success */
	return 0;
}


error_t write(char *file_name, int N_conf, double *MSD)
{
	printf("Writing output to '%s'...\n", file_name);


	/* Opening the output file */
	FILE* output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening the output file");
		return EIO;
	}


	/* Actually writing */
	// Writing the header
	fprintf(output, "# tau\tMSD\n");

	for (int tau = 0 ; tau < N_conf ; tau++)
		fprintf(output, " %d\t%lf\n", tau, MSD[tau]);
	
	
	/* Success */
	// Closing the output file
	fclose(output);

	// Exiting normally
	return 0;
}


int main(int argc, char **argv)
{
	/* Parsing the arguments */
	Arguments arguments;

	// Options' default values
	set_default_options(&arguments);

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
	{
		perror("Parsing");
		exit(EXIT_FAILURE);
	}


	/* Reading the configurations */
	int N_conf, *N_selection, *steps;
	double **bounds;
	Atom **atoms;

	// Reading the file
	if ((errno = read_trajectory(&arguments, &N_conf, &steps, &N_selection, &bounds, &atoms)) != 0)
		goto READ;

	
	/* Computing the MSD */
	double *MSD;
	if ((errno = compute_msd(N_conf, N_selection, atoms, &MSD)) != 0)
		goto READ;


	/* Writing the full MSD */
	if ((errno = write("output/msd.dat", N_conf, MSD)) != 0)
		goto COMPUTE;


	/* Success */
	exit(EXIT_SUCCESS);


	/* Errors */
	COMPUTE: free(MSD);
	READ: free(atoms), free(bounds), free(steps), free(N_selection);
	exit(EXIT_FAILURE);
}
