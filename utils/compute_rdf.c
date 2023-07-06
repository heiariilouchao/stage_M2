# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <argp.h>
# include <math.h>
# define STR_BUFF_LIMIT 256
# define STR_LABELS_LIMIT 8


static char doc[] = "Computing the RDFs and MSDs of a .lammpstrj configurations file.";

static char args_doc[] = "CONF_FILE";

static struct argp_option options[] =
{
    {
        "bins",
        'b',
        "BINS",
        OPTION_ARG_OPTIONAL,
        "The number of bins to compute the RDFs. Defaults to 50."
    },
    {
        "steps",
        's',
        "STEPS",
        0,
        "The maximum number of configurations to process."
    },
    {
        "skip",
        'k',
        "SKIP",
        OPTION_ARG_OPTIONAL,
        "The number of steps to skip for the RDFs. Defaults to 0."
    },
    {0}
};

struct arguments
{
    char *args[1];
    int N_bins, N_maxsteps, skip;
};

static error_t parse(int key, char *arg, struct argp_state *state)
{
    struct arguments *args = state->input;
    
    switch (key)
    {
        case 'b':
            args->N_bins = atoi(arg);
            if (args->N_bins <= 0)
                return EINVAL;
            break;
        case 's':
            args->N_maxsteps = atoi(arg);
            if (args->N_maxsteps <= 0)
                return EINVAL;
            break;
        case 'k':
            args->skip = atoi(arg);
            if (args->skip < 0)
                return EINVAL;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= 1)
                argp_usage(state);
            args->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num <1)
                argp_usage(state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp parser =
{
    options,
    parse,
    args_doc,
    doc
};


error_t read_parameters(char* file_name, int* N_conf, int* N_atoms, double bounds[3][2], char ***labels, int *N_labels)
{
    FILE* conf_file_ptr = fopen(file_name, "r");

    if (conf_file_ptr == NULL)
    {
        perror("Opening the configuration file to read the parameters");
        exit(EXIT_FAILURE);
    }

    char dump[STR_BUFF_LIMIT], str_value[STR_BUFF_LIMIT];
    int error, int_value, count = 0;
    bool atoms_verif = true, bounds_verif = true;

    *labels = (char **) malloc(sizeof(char **));
    if (*labels == NULL)
    {
        fclose(conf_file_ptr);
        perror("Allocating the labels");
        return ENOMEM;
    }

    int n = 0;

    while (fgetc(conf_file_ptr) != EOF)
    {
        ungetc('I', conf_file_ptr);
        error = fscanf(conf_file_ptr, "ITEM: %s", str_value);
        
        if (strcmp(str_value, "NUMBER") == 0 && atoms_verif)
        {
            if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
            {
                fclose(conf_file_ptr);
                perror("Dumping the rest of the line (NUMBER)");
                return EIO;
            }

            if (fscanf(conf_file_ptr, "%d", &int_value) != 1)
            {
                fclose(conf_file_ptr);
                perror("Number of atoms");
                return EIO;
            }

            *N_atoms = int_value;
            atoms_verif = false;
        }
        else if (strcmp(str_value, "BOX") == 0 && bounds_verif)
        {
            if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
            {
                fclose(conf_file_ptr);
                perror("Dumping the rest of the line (BOX)");
                return EIO;
            }

            for (int i = 0 ; i < 3 ; i++)
            {
                    if (fscanf(conf_file_ptr, "%lf %lf", &(bounds[i][0]), &(bounds[i][1])) != 2)
                    {
                        fclose(conf_file_ptr);
                        perror("Box boundaries");
                        return EIO;
                    }
            }
            bounds_verif = false;
        }
        else if (strcmp(str_value, "ATOMS") == 0)
        {
            count++;

            if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
            {
                fclose(conf_file_ptr);
                perror("Dumping the rest of the line (ATOMS)");
                return EIO;
            }

            for (int i = 0 ; i < *N_atoms ; i++)
            {
                if (fscanf(conf_file_ptr, "%*d %s", str_value) != 1)
                {
                    fclose(conf_file_ptr);
                    perror("Atom labels");
                    return EIO;
                }

                if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
                {
                    fclose(conf_file_ptr);
                    perror("Dumping the rest of the line (ATOMS LABELS)");
                    return EIO;
                }

                if (n == 0)
                {
                    (*labels)[n] = (char *) malloc(strlen(str_value) + 1);

                    if ((*labels)[n] == NULL)
                    {
                        fclose(conf_file_ptr);
                        perror("Labels first slot");
                        return ECANCELED;
                    }

                    strcpy((*labels)[n], str_value);
                    n++;
                }
                else
                {
                    bool verif = false;
                    for (int j = 0 ; j < n ; j++)
                    {
                        verif = (strcmp(str_value, (*labels)[j]) == 0);
                        if (verif)
                            break;
                    }
                    if (!verif)
                    {
                        (*labels)[n] = (char *) malloc(strlen(str_value) + 1);

                        if ((*labels)[n] == NULL)
                        {
                            fclose(conf_file_ptr);
                            perror("Labels later slot");
                            return ECANCELED;
                        }

                        strcpy((*labels)[n], str_value);
                        n++;
                    }
                }
            }
        }
    }

    *N_conf = count;
    *N_labels = n;
    fclose(conf_file_ptr);
    return 0;
}


error_t compute_pairs(int N_labels, char **labels, char ***pairs)
{
    *pairs = (char **) malloc(sizeof(char **));

    if (*pairs == NULL)
    {
        perror("Pairs first attempt");
        return ENOMEM;
    }

    for (int l1 = 0, p = 0 ; l1 < N_labels ; l1++)
        for (int l2 = l1 ; l2 < N_labels ; l2++, p++)
        {
            (*pairs)[p] = (char *) malloc(strlen(labels[l1]) + strlen(labels[l2]) + 1);

            if ((*pairs)[p] == NULL)
            {
                perror("Pairs on later attempt");
                return ECANCELED;
            }

            strcpy((*pairs)[p], labels[l1]);
            strcat((*pairs)[p], labels[l2]);
        }

    return 0;
}


int main(int argc, char **argv)
{
    /* Error code */
    error_t err;


    /* Parsing the arguments */
    struct arguments arguments;

    // Options' default values
    arguments.N_bins = 50;
    arguments.N_maxsteps = 50;
    arguments.skip = 0;

    // Actually parsing
    err = argp_parse(&parser, argc, argv, 0, 0, &arguments);

    // Error handling
    if (err != 0)
    {
        perror("Parsing");
        goto ERROR;
    }


    /* Reading the parameters */
    int N_conf, N_atoms, N_labels;
    double bounds[3][2];
    char **labels;

    // Actually reading
    err = read_parameters(arguments.args[0], &N_conf, &N_atoms, bounds, &labels, &N_labels);

    // Error handling
    if (err != 0)
    {
        if (err == ENOMEM)
            goto ERROR;
        else if (err == EIO || err == ECANCELED)
            goto LABELS;
    }


    /* Computing extra parameters */
    double hbounds[3][2], cutoff, delta, V = 1.;
    for (int d = 0 ; d < 3 ; d++)
    {
        hbounds[d][0] = bounds[d][0] / 2.;
        hbounds[d][1] = bounds[d][1] / 2.;

        V *= bounds[d][1] - bounds[d][0];

        if (d == 0)
            cutoff = hbounds[d][1] - hbounds[d][0];
        else
            cutoff = (hbounds[d][1] - hbounds[d][0] < cutoff) ? hbounds[d][1] - hbounds[d][0] : cutoff;
    }
    delta = cutoff / arguments.N_bins;


    /* Computing the pairs */
    int N_pairs = N_labels * (N_labels + 1) / 2;
    char **pairs;
    err = compute_pairs(N_labels, labels, &pairs);

    // Error handling
    if (err != 0)
    {
        if (err == ENOMEM)
            goto LABELS;
        else if (err == ECANCELED)
            goto PAIRS;
    }


    /* Declaring and allocating the arrays */
    // The pairs-dependent
    int **histograms = (int **) malloc(N_pairs * sizeof(int *));

    if (histograms == NULL)
    {
        perror("Allocating the histograms");
        goto PAIRS;
    }

    double **RDFs = (double **) malloc(N_pairs * sizeof(double *));

    if (RDFs == NULL)
    {
        perror("Allocating the RDFs");
        goto HIST;
    }

    for (int p = 0 ; p < N_pairs ; p++)
    {
        histograms[p] = (int *) calloc(arguments.N_bins, sizeof(int));
        if (histograms[p] == NULL)
        {
            perror("Slot of histogram");
            goto RDF;
        }
        
        RDFs[p] = (double *) malloc(arguments.N_bins * sizeof(double));
        if (RDFs[p] == NULL)
        {
            perror("Slot of RDF");
            goto RDF;
        }
    }


    // The atoms-dependent
    double **current_positions = (double **) malloc(N_atoms * sizeof(double *));
    if (current_positions == NULL)
    {
        perror("Allocating the current positions");
        goto RDF;
    }

    for (int a = 0 ; a < N_atoms ; a++)
    {
        current_positions[a] = (double *) malloc(3 * sizeof(double));
        if (current_positions[a] == NULL)
        {
            perror("Slot of current position");
            goto CPOSITIONS;
        }
    }


    /* Building the hash-tables */
    // Declaring and allocating the arrays
    int **indices = (int **) malloc(N_labels * sizeof(int *));
    if (indices == NULL)
    {
        perror("Allocating the indices");
        goto CPOSITIONS;
    }

    int *n_labels = (int *) calloc(N_labels, sizeof(int));
    if (n_labels == NULL)
    {
        perror("Allocating the number of labels");
        goto INDICES;
    }

    for (int l = 0 ; l < N_labels ; l++)
    {
        indices[l] = (int *) malloc(N_atoms * sizeof(int));
        if (indices[l] == NULL)
        {
            perror("Slot of index");
            goto NLABELS;
        }
    }

    // Opening the configuration file
    FILE* conf_file_ptr = fopen(arguments.args[0], "r");
    char dump[STR_BUFF_LIMIT];

    // Error handling
    if (conf_file_ptr == NULL)
    {
        perror("Reading the configuration file to build the hash-tables");
        goto NLABELS;
    }

    // Dumping the first lines
    err = 0;
    for (int l = 0 ; l < 9 ; l++)
        err += (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);
    
    if (err != 0)
    {
        perror("Dumping the first lines of the initial configuration");
        goto CONF_FILE;
    }

    // Reading the first configuration
    for (int a = 0 ; a < N_atoms ; a++)
    {
        // Declaring the temporary variables
        int atom_index;
        char label[STR_LABELS_LIMIT];

        // Reading the first two fields of a line
        if (fscanf(conf_file_ptr, "%d %s", &atom_index, label) != 2)
        {
            errno = EIO;
            perror("Reading the first fields of a line");
            goto CONF_FILE;
        }

        if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
        {
            errno = EIO;
            perror("Dumping the rest of the line");
            goto CONF_FILE;
        }

        atom_index--;
        for (int l = 0 ; l < N_labels ; l++)
            if (strcmp(labels[l], label) == 0)
            {
                indices[l][n_labels[l]] = atom_index;
                n_labels[l]++;
            }
    }

    // Resizing the arrays
    for (int l = 0 ; l < N_labels ; l++)
    {
        indices[l] = (int *) realloc(indices[l], n_labels[l] * sizeof(int));
        if (indices[l] == NULL)
        {
            perror("Resizing the indices arrays");
            goto CONF_FILE;
        }
    }


    /* Processing the configurations */
    // Opening the configuration file
    conf_file_ptr = freopen(arguments.args[0], "r", conf_file_ptr);

    // Error handling
    if (conf_file_ptr == NULL)
    {
        perror("Reading the configuration file to process it");
        goto NLABELS;
    }

    // Going through the file
    for (int step = 0 ; step < arguments.N_maxsteps ; step++)
    {
        if (step % 100 == 0)
            printf("\rstep: %d / %d", step + 1, arguments.N_maxsteps);
        // Declaring the temporary variables
        int atom_index;

        // Dumping the first lines
        for (int l = 0 ; l < 9 ; l++)
            if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
            {
                err = EIO;
                perror("Dumping the first lines of the current configuration");
                goto CONF_FILE;
            }
        
        // Reading the current configuration
        for (int a = 0 ; a < N_atoms ; a++)
        {
            if (fscanf(conf_file_ptr, "%d %*s", &atom_index) != 1)
            {
                err = EIO;
                perror("Reading the atom index");
                goto CONF_FILE;
            }

            atom_index--;
            if (fscanf(conf_file_ptr, "%lf %lf %lf", &(current_positions[atom_index][0]),
                                                     &(current_positions[atom_index][1]),
                                                     &(current_positions[atom_index][2])) != 3)
            {
                err = EIO;
                perror("Reading the positions");
                goto CONF_FILE;
            }
        
            if (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL)
            {
                err = EIO;
                perror("Dumping the rest of the line");
                goto CONF_FILE;
            }
        }

        if (step >= arguments.skip)
        {
            // Incrementing the histograms
            for (int l1 = 0, h = 0 ; l1 < N_labels ; l1++)
            {
                for (int l2 = l1 ; l2 < N_labels ; l2++, h++)
                {
                    if (l1 != l2)
                        for (int i1 = 0 ; i1 < n_labels[l1] ; i1++)
                        {
                            for (int i2 = 0 ; i2 < n_labels[l2] ; i2++)
                            {
                                double distance = 0.;
                                for (int d = 0 ; d < 3 ; d++)
                                {
                                    double diff = current_positions[indices[l2][i2]][d] - current_positions[indices[l1][i1]][d];
                                    if (diff < -hbounds[d][1])
                                        diff += bounds[d][1];
                                    else if (diff > hbounds[d][1])
                                        diff -= bounds[d][1];
                                    distance += diff * diff;
                                }
                                distance = sqrt(distance);
                                int index = (int) (distance / delta);
                                if (index < arguments.N_bins)
                                    histograms[h][index]++;
                            }
                        }
                    else
                        for (int i1 = 0 ; i1 < n_labels[l1] ; i1++)
                        {
                            for (int i2 = i1 + 1 ; i2 < n_labels[l1] ; i2++)
                            {
                                double distance = 0.;
                                for (int d = 0 ; d < 3 ; d++)
                                {
                                    double diff = current_positions[indices[l1][i2]][d] - current_positions[indices[l1][i1]][d];
                                    if (diff < -hbounds[d][1])
                                        diff += bounds[d][1];
                                    else if (diff > hbounds[d][1])
                                        diff -= bounds[d][1];
                                    distance += diff * diff;
                                }
                                distance = sqrt(distance);
                                int index = (int) (distance / delta);
                                if (index < arguments.N_bins)
                                    histograms[h][index] += 2;
                            }
                        }
                }
            }
        }
    }


    /* Freeing the pointers */
    fclose(conf_file_ptr);
    free(indices);
    free(current_positions);


    /* Averaging the RDFs */
    // Computing the distance range
    double *r = (double *) malloc(arguments.N_bins * sizeof(double));
    if (r == NULL)
    {
        perror("Allocating the distance range array");
        goto NLABELS2;
    }

    for (int b = 0 ; b < arguments.N_bins ; b++)
        r[b] = (double) b * delta;

    // Computing the number of atoms for each pair
    int *N_atoms_pairs = (int *) malloc(N_pairs * sizeof(int));
    if (N_atoms_pairs == NULL)
    {
        perror("Allocating the numbers of atoms per pairs array");
        goto R;
    }

    for (int l1 = 0, p = 0 ; l1 < N_labels ; l1++)
        for (int l2 = l1 ; l2 < N_labels ; l2++, p++)
            N_atoms_pairs[p] = (l1 != l2) ? n_labels[l1] * n_labels[l2] : n_labels[l1] * (n_labels[l1] - 1);

    // Normalizing and averaging the histograms
    const double constant = 1. / (4. / 3. * M_PI / V * (arguments.N_maxsteps - arguments.skip));
    for (int p = 0 ; p < N_pairs ; p++)
        for (int b = 0 ; b < arguments.N_bins ; b++)
            RDFs[p][b] = constant * histograms[p][b] / (N_atoms_pairs[p] * (pow(r[b] + delta, 3) - pow(r[b], 3)));


    /* Freeing the unused pointers */
    free(n_labels);
    free(N_atoms_pairs);


    /* Shifting the distance range */
    for (int b = 0 ; b < arguments.N_bins ; b++)
        r[b] += delta / 2.;


    /* Writing the results */
    FILE *rdfs_file_ptr;

    // Writing the RDFs
    rdfs_file_ptr = fopen("rdf.dat", "w");

    // Error handling
    if (rdfs_file_ptr == NULL)
    {
        perror("Writing the RDFs");
        goto R2;
    }

    fprintf(rdfs_file_ptr, "# %7s", "r");
    
    for (int p = 0 ; p < N_pairs ; p++)
    {
        char *str = (char *) malloc(strlen(pairs[p]) + 3);
        strcpy(str, "g_");
        fprintf(rdfs_file_ptr, " %7s", strcat(str, pairs[p]));
    }
    fprintf(rdfs_file_ptr, "\n");

    for (int b = 0 ; b < arguments.N_bins ; b++)
    {
        fprintf(rdfs_file_ptr, "  %2.4lf", r[b]);
        for (int p = 0 ; p < N_pairs ; p++)
            fprintf(rdfs_file_ptr, " %2.4lf", RDFs[p][b]);
        fprintf(rdfs_file_ptr, "\n");
    }
    
    fclose(rdfs_file_ptr);


    /* Exiting normally */
    free(r);
    free(RDFs);
    free(histograms);
    free(pairs);
    free(labels);
    exit(EXIT_SUCCESS);


    /* Freeing the pointers and exiting */
    CONF_FILE: fclose(conf_file_ptr);
    NLABELS: free(n_labels);
    INDICES: free(indices);
    CPOSITIONS: free(current_positions);
    RDF: free(RDFs);
    HIST: free(histograms);
    PAIRS: free(pairs);
    LABELS: free(labels);
    ERROR: errno = err ;
    exit(EXIT_FAILURE);


    /* Freeing the second batch of pointers and exiting */
    NATOMSPAIRS: free(N_atoms_pairs);
    R: free(r);
    NLABELS2: free(n_labels);
    free(RDFs);
    free(histograms);
    free(pairs);
    free(labels);
    errno = err;
    exit(EXIT_FAILURE);


    /* Freeing the third batch of pointers and exiting */
    R2: free(r);
    free(RDFs);
    free(histograms);
    free(pairs);
    free(labels);
    errno = err;
    exit(EXIT_FAILURE);
}