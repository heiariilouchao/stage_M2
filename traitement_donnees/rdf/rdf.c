# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <math.h>

# include "../utils.h"
# include "../parse/parse.h"
# include "../read/read.h"


char doc[] = "Computing the RDFs of a .lammpstrj configurations file.";

char args_doc[] = "CONF_FILE";


int compute_pairs(int N_elements, int ***map)
{
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


double compute_cutoff(double **bounds)
{
	printf("Computing the cutoff...\n");


	double cutoff = bounds[0][1] - bounds[0][0];
	if (bounds[0][3] - bounds[0][2] < cutoff)
		cutoff = bounds[0][3] - bounds[0][2];
	if (bounds[0][5] - bounds[0][4] < cutoff)
		cutoff = bounds[0][5] - bounds[0][4];
	

	return floor(cutoff / 2.);
}


int compute_rdf(int N_conf, Arguments arguments, int *N_selection, double **bounds, Atom **atoms, double **r, double ***RDF)
{
	int N_elements = arguments.N_elements;
	int N_bins = arguments.N_bins;
	int N_pairs = arguments.N_pairs;
	double cutoff = compute_cutoff(bounds);
	double delta = cutoff / N_bins;

	printf("Cutoff informations:\n\tcutoff: %lf\n\tbins: %d\n\tdelta: %lf\n", cutoff, N_bins, delta);

	printf("Pairs informations:\n\tN: %d\n\tpairs:", N_pairs);
	for (int e1 = 0, p = 0 ; e1 < N_elements ; e1++)
		for (int e2 = e1 ; e2 < N_elements ; e2++, p++)
			printf(" %s%s", arguments.elements[e1], arguments.elements[e2]);
	printf("\n");


	/* Allocating the arrays */
	int *N, **map, **hist;

	if (((*r) = malloc(N_bins * sizeof(double))) == NULL)
	{
		perror("Allocating an array (r)");
		return ENOMEM;
	}

	if ((N = calloc(N_pairs, sizeof(int))) == NULL)
	{
		perror("Allocating an array (N)");
		goto R;
	}

	if (compute_pairs(N_elements, &map) != 0)
		goto N;

	if (((hist = malloc(N_pairs * sizeof(int *)))) == NULL)
	{
		perror("Allocating an array (hist)");
		goto MAP;
	}

	if (((*RDF) = malloc(N_pairs * sizeof(double *))) == NULL)
	{
		perror("Allocating an array (RDF)");
		goto HIST;
	}

	for (int p = 0 ; p < N_pairs ; p++)
	{
		if ((hist[p] = calloc(N_bins, sizeof(int))) == NULL)
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
	{
		printf("conf: %d / %d\r", c + 1, N_conf);
		for (int i = 0 ; i < N_selection[c] ; i++)
			for (int j = i + 1 ; j < N_selection[c] ; j++)
			{
				int p = map[atoms[c][i].element_ID][atoms[c][j].element_ID];
				if (c == 0)
					N[p] += 2;

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
					hist[p][bin] += 2;
			}
	}
	

	/* Computing the RDFs */
	printf("Computing the RDFs...\n");

	double V = 1.;
	for (int d = 0 ; d < 3 ; d++)
		V *= bounds[0][2 * d + 1] - bounds[0][2 * d];
	double coeff = 1. / (4. / 3. * M_PI / V);
	
	for (int p = 0 ; p < N_pairs ; p++)
		for (int b = 0 ; b < N_bins ; b++)
			(*RDF)[p][b] = (double) coeff * hist[p][b] / N_conf / (N[p] * (pow((*r)[b] + delta, 3) - pow((*r)[b], 3)));

	// Shifting the radial range
	for (int b = 0 ; b < N_bins ; b++)
		(*r)[b] += delta / 2.;
	

	/* Success */
	free(hist);
	free(map);
	free(N);

	// Exiting normally
	return 0;


	/* Errors */
	RDF: free(RDF);
	HIST: free(hist);
	MAP: free(map);
	N: free(N);
	R: free(r);
	return ENOMEM;
}


int write(char *file_name, Arguments arguments, double *r, double **RDF)
{
	printf("Writing the output to '%s'...\n", file_name);

	int N_elements = arguments.N_elements;
	char **elements = arguments.elements;
	int N_bins = arguments.N_bins;
	int N_pairs = arguments.N_pairs;


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
	Arguments arguments;

	// Options' default values
	set_default_options(&arguments);

	// Actually parsing
	if (argp_parse(&parser, argc, argv, 0, 0, &arguments) != 0)
		exit(EXIT_FAILURE);


	/* Reading the configurations */
	int N_conf, *N_selection, *steps;
	double **bounds;
	Atom **atoms;

	// Reading the file
	if ((errno = read_trajectory(&arguments, &N_conf, &steps, &N_selection, &bounds, &atoms)) != 0)
		exit(EXIT_FAILURE);
	

	/* Computing the RDFs */
	double *r, **RDF;
	if ((errno = compute_rdf(N_conf, arguments, N_selection, bounds, atoms, &r, &RDF)) != 0)
		goto READ;
	

	/* Writing the output */
	if ((errno = write("output/rdf.dat", arguments, r, RDF)) != 0)
		goto RDF;


	/* Success */
	exit(EXIT_SUCCESS);


	/* Error handling */
	RDF: free(RDF), free(r);
	READ: free(atoms), free(bounds), free(steps), free(N_selection);
	exit(EXIT_FAILURE);
}
