# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <math.h>

# include "../utils/utils.h"
# include "../parse/parse.h"
# include "../read/read.h"


int compute_pairs(int N_elements, char **elements, int *N_pairs, char ***pairs)
{
	printf("Computing the pairs...\n");
	

	/* Allocating the array */
	*N_pairs = N_elements * (N_elements + 1) / 2;
	if ((*pairs = malloc(*N_pairs * sizeof(char *))) == NULL)
	{
		perror("Allocating an array (pairs)");
		return ENOMEM;
	}

	for (int p = 0 ; p < *N_pairs ; p++)
		if (((*pairs)[p] = malloc((2 * STR_ELEMENT_LIMIT + 2) * sizeof(char))) == NULL)
		{
			perror("Allocating an array slot (pairs[])");
			goto NOMEM;
		}


	/* Computing the pairs */
	for (int e1 = 0, p = 0 ; e1 < N_elements ; e1++)
		for (int e2 = e1 ; e2 < N_elements ; e2++, p++)
		{
			char pair[2 * STR_ELEMENT_LIMIT + 2] = {0};
			strcpy(pair, elements[e1]);
			strcat(pair, ",");
			strcat(pair, elements[e2]);
			strcpy((*pairs)[p], pair);
		}
	
	/* Prompting informations */
	printf("Pairs informations:\n\tN_pairs: %d\n\tpairs:", *N_pairs);
	for (int p = 0 ; p < *N_pairs ; p++)
		printf(" %s", (*pairs)[p]);
	printf("\n");


	/* Success */
	// Exiting normally
	return 0;

	
	/* Errors */
	NOMEM: free(pairs);
	return ENOMEM;
}


double compute_cutoff(Box *box)
{
	printf("Computing the cutoff...\n");


	double length = box[0].x_max - box[0].x_min;
	if (box[0].y_max - box[0].y_min < length)
		length = box[0].y_max - box[0].y_min;
	if (box[0].z_max - box[0].z_min < length)
		length = box[0].z_max - box[0].z_min;

	return floor(length / 2.);
}


int compute_rdf(int N_configurations, Box *box, int N_bins, bool are_identical, int *N1, Atom **a1, int *N2, Atom **a2, double **r, double **rdf)
{
	printf("Computing the RDF...\n");


	/* Computing the parameters */
	double cutoff = compute_cutoff(box);
	double delta = (double) cutoff / N_bins;


	/* Allocating the arrays */
	int *hist;

	if ((*r = malloc(N_bins * sizeof(double))) == NULL)
	{
		perror("Allocating an array (r)");
		return ENOMEM;
	}

	if ((hist = calloc(N_bins, sizeof(int))) == NULL)
	{
		perror("Allocating an array (hist)");
		goto R;
	}

	if ((*rdf = malloc(N_bins * sizeof(double))) == NULL)
	{
		perror("Allocating an array (RDF)");
		goto HIST;
	}

	for (int b = 0 ; b < N_bins ; b++)
		(*r)[b] = delta * b;
	

	/* Incrementing the histogram */
	printf("Incrementing the histogram...\n");
	for (int c = 0 ; c < N_configurations ; c++)
	{
		printf("conf: %d / %d\r", c, N_configurations);
		for (int i = 0 ; i < N1[c] ; i++)
			for (int j = 0 ; j < N2[c] ; j++)
			{
				if (are_identical && a1[c][i].serial == a2[c][j].serial)
					continue;
				double r2 = 0.;

				double length, diff;
				length = box[c].x_max - box[c].x_min;
				diff = a2[c][j].x - a1[c][i].x;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				length = box[c].y_max - box[c].y_min;
				diff = a2[c][j].y - a1[c][i].y;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;
				
				length = box[c].z_max - box[c].z_min;
				diff = a2[c][j].z - a1[c][i].z;
				if (diff < - length / 2.)
					diff += length;
				else if (length / 2. < diff)
					diff -= length;
				r2 += diff * diff;

				int bin = (int) (sqrt(r2) / delta);
				if (bin < N_bins)
					hist[bin]++;
			}
	}


	/* Computing the RDF */
	// Prompting
	printf("Computing the RDF...\n");

	// Computing the parameters
	double V = (box[0].x_max - box[0].x_min) * (box[0].y_max - box[0].y_min) * (box[0].z_max - box[0].z_min);
	
	int N = N1[0];
	if (are_identical)
		N *= N1[0] - 1;
	else
		N *= N2[0];
	
	double coeff = 1. / (4. / 3. * M_PI * N / V);
	
	// Actually computing the RDF
	for (int b = 0 ; b < N_bins ; b++)
		(*rdf)[b] = (double) coeff * hist[b] / N_configurations / (pow((*r)[b] + delta, 3) - pow((*r)[b], 3));

	// Shifting the radial range
	for (int b = 0 ; b < N_bins ; b++)
		(*r)[b] += delta / 2.;
	

	/* Success */
	free(hist);

	// Exiting normally
	return 0;
	

	/* Errors */
	HIST: free(hist);
	R: free(r);
	return ENOMEM;
}


int write_pairs(char *file_name, int N_pairs, char **pairs, int N_bins, double *r, double **rdf)
{
	printf("Writing the output to '%s'...\n", file_name);


	/* Opening the file */
	FILE* output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing */
	fprintf(output, "# r");
	for (int p = 0 ; p < N_pairs ; p++)
		fprintf(output, " g_{%s}", pairs[p]);
	fprintf(output, "\n");


	for (int b = 0 ; b < N_bins ; b++)
	{
		fprintf(output, "  %lf", r[b]);
		for (int p = 0 ; p < N_pairs ; p++)
			fprintf(output, " %lf", rdf[p][b]);
		fprintf(output, "\n");
	}


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}


int write_rdf(char *file_name, char *label, int N_bins, double *r, double *rdf)
{
	printf("Writing the output to '%s'...\n", file_name);


	/* Opening the file */
	FILE* output;
	if ((output = fopen(file_name, "w")) == NULL)
	{
		perror("Opening a file (output)");
		return EIO;
	}


	/* Writing */
	fprintf(output, "# r g_{%s}\n", label);

	for (int b = 0 ; b < N_bins ; b++)
		fprintf(output, "  %lf %lf\n", r[b], rdf[b]);


	/* Success */
	// Closing the file
	fclose(output);

	// Exiting normally
	return 0;
}
