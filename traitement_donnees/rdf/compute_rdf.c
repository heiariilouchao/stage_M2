# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <math.h>

# include "utils.h"
# include "parse.h"
# include "read.h"


int compute_pairs(int N_elements, int ***map)
{
	printf("Computing the pairs...\n");


	/* Allocating the arrays */
	if ((*map = malloc(N_elements * sizeof(int *))) == NULL)
	{
		perror("Allocating an array (map)");
		return ENOMEM;
	}

	for (int e = 0 ; e < N_elements ; e++)
		if (((*map)[e] = calloc(N_elements, sizeof(int))) == NULL)
		{
			perror("Allocating an array slot (map[])");
			goto MAP;
		}
	

	/* Computing the pairs */
	for (int i = 0, p = 0 ; i < N_elements ; i++)
		for (int j = i ; j < N_elements ; j++, p++)
		{
			(*map)[i][j] = p;
			(*map)[j][i] = p;
		}


	/* Success */
	// Exiting normally
	return 0;


	/* Errors */
	MAP: free(map);
	return ENOMEM;
}


int compute_rdf(int N_conf, int N_elements, int *N_selection, double **bounds, int N_bins, double cutoff, int N_pairs, int **map, Atom **atoms, double **r, double ***RDF)
{
	double delta = cutoff / N_bins;

	/* Allocating the arrays */
	int **hist;

	if (((*r) = malloc(N_bins * sizeof(double))) == NULL)
	{
		perror("Allocating an array (r)");
		return ENOMEM;
	}

	if (((hist = malloc(N_pairs * sizeof(double *)))) == NULL)
	{
		perror("Allocating an array (hist)");
		goto R;
	}

	if (((*RDF) = malloc(N_pairs * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (RDF)");
		goto HIST;
	}

	for (int p = 0 ; p < N_pairs ; p++)
	{
		if ((hist[p] = calloc(N_bins, sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (hist[])");
			goto RDF;
		}

		if (((*RDF)[p] = calloc(N_bins, sizeof(double))) == NULL)
		{
			perror("Allocating an array slot (RDF[])");
			goto RDF;
		}
	}

	// Initializing the distance range
	for (int b = 0 ; b < N_bins ; b++)
		(*r)[b] = b * delta;


	/* Incrementing the histograms */
	printf("Incrementing the histograms...\n");

	for (int c = 0 ; c < N_conf ; c++)
		for (int i = 0 ; i < N_selection[c] ; i++)
			for (int j = i + 1 ; j < N_selection[c] ; j++)
			{
				int p = map[atoms[c][i].element_ID][atoms[c][j].element_ID];

				double r2 = 0.;

				double length, diff;
				length = bounds[c][1] - bounds[c][0];
				diff = atoms[c][j].x - atoms[c][i].x;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				length = bounds[c][3] - bounds[c][2];
				diff = atoms[c][j].y - atoms[c][i].y;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;
				
				length = bounds[c][5] - bounds[c][4];
				diff = atoms[c][j].z - atoms[c][i].z;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				int bin = (int) (sqrt(r2) / delta);
				if (bin < N_bins)
					hist[p][bin] += (atoms[c][i].element_ID == atoms[c][j].element_ID) ? 2 : 1;
			}
	

	/* Computing the RDFs */
	printf("Normalizing the RDFs...\n");

	double V = 1.;
	for (int d = 0 ; d < 3 ; d++)
		V *= bounds[0][2 *d + 1] - bounds[0][2 * d];
	double coeff = 1. / (4. / 3. * M_PI / V);

	for (int e1 = 0, p = 0 ; e1 < N_elements ; e1++)
	{
		for (int e2 = e1 ; e2 < N_elements ; e2++, p++)
		{
			int N1 = 0, N2 = 0, N;
			if (e1 == e2)
			{
				for (int a = 0 ; a < N_selection[0] ; a++)
					if (atoms[0][a].element_ID == e1)
						N1++;
				N = N1 * (N1 - 1);
			}
			else
			{
				for (int a = 0 ; a < N_selection[0] ; a++)
					if (atoms[0][a].element_ID == e1)
						N1++;
					else if (atoms[0][a].element_ID == e2)
						N2++;
				N = N1 * N2;
			}

			for (int b = 0 ; b < N_bins ; b++)
				(*RDF)[p][b] = (double) coeff * hist[p][b] / N_conf / (N * (pow((*r)[b] + delta, 3) - pow((*r)[b], 3)));
		}
	}
	

	/* Success */
	return 0;


	/* Errors */
	RDF: free(RDF);
	HIST: free(hist);
	R: free(r);
	return ENOMEM;
}


int write(char *file_name, int N_elements, char **elements, int N_bins, int N_pairs, double *r, double **RDF)
{
	printf("Writing the output...\n");

	/* Opening the file */
	FILE* output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing */
	fprintf(output, "# r");
	for (int e1 = 0 ; e1 < N_elements ; e1++)
		for (int e2 = e1 ; e2 < N_elements ; e2++)
			fprintf(output, " g_%s%s", elements[e1], elements[e2]);
	fprintf(output, "\n");


	for (int b = 0 ; b < N_bins ; b++)
	{
		fprintf(output, "  %lf", r[b]);
		for (int p = 0 ; p < N_pairs ; p++)
			fprintf(output, " %lf", RDF[p][b]);
		fprintf(output, "\n");
	}


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}


int main(int argc, char **argv)
{
	/* Parsing the arguments */
	struct arguments arguments;

	// Options' default values
	arguments.start = 0;
	arguments.N_bins = 100;

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
		exit(EXIT_FAILURE);


	/* Reading the configurations */
	// Declaring the arrays
	int N_conf, *N_selection, *steps;
	double **bounds;
	Atom **atoms;

	// Reading the file
	if ((errno = read_trajectory(arguments.args[0], arguments.start, arguments.N_elements, arguments.labels, arguments.elements, &N_conf, &steps, &N_selection, &bounds, &atoms)) != 0)
		exit(EXIT_FAILURE);
	

	/* Computing the pairs */
	int **pair_map;
	if ((errno = compute_pairs(arguments.N_elements, &pair_map)) != 0)
		goto READ;
	
	for (int i = 0 ; i < arguments.N_elements ; i++)
		for (int j = 0 ; j < arguments.N_elements ; j++)
			printf("%d %d %d\n", i, j, pair_map[i][j]);
	printf("\n");
	

	/* Computing the RDFs */
	double *r, **RDF;
	double cutoff = 10.0;
	if ((errno = compute_rdf(N_conf, arguments.N_elements, N_selection, bounds, arguments.N_bins, cutoff, arguments.N_pairs, pair_map, atoms, &r, &RDF)) != 0)
		goto PAIRS;
	

	/* Writing the output */
	if ((errno = write("output/rdf.dat", arguments.N_elements, arguments.elements, arguments.N_bins, arguments.N_pairs, r, RDF)) != 0)
		goto RDF;


	/* Exiting normally */
	exit(EXIT_SUCCESS);


	/* Error handling */
	RDF: free(RDF), free(r);
	PAIRS: free(pair_map);
	READ: free(atoms), free(bounds), free(steps), free(N_selection);
	exit(EXIT_FAILURE);
}
