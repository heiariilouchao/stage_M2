# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <string.h>
# include <errno.h>
# include <argp.h>
# include <math.h>
# include <omp.h>
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
        0,
        "The number of bins to compute the RDFs."
    },
    {
        "steps",
        's',
        "STEPS",
        0,
        "The maximum number of configurations to process."
    },
    {0}
};

struct arguments
{
    char *args[1];
    int N_bins, N_maxsteps;
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
        perror("Configuration file");
        exit(EXIT_FAILURE);
    }

    char dump[STR_BUFF_LIMIT], str_value[STR_BUFF_LIMIT];
    int error, int_value, count = 0;
    double double_value1, double_value2;
    bool atoms_verif = true, bounds_verif = true, char_verif = true;
    *labels = (char **) malloc(sizeof(char **));
    if (*labels == NULL)
    {
        fclose(conf_file_ptr);
        perror("Labels");
        return ENOMEM;
    }

    int n = 0;

    while (fgetc(conf_file_ptr) != EOF)
    {
        ungetc('I', conf_file_ptr);
        error = fscanf(conf_file_ptr, "ITEM: %s", str_value);
        
        if (strcmp(str_value, "NUMBER") == 0 && atoms_verif)
        {
            error = (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);
            error = fscanf(conf_file_ptr, "%d", &int_value);
            if (error == 1)
            {
                *N_atoms = int_value;
                atoms_verif = false;
            }
            else
            {
                fclose(conf_file_ptr);
                perror("Number of atoms");
                return EIO;
            }
        }
        else if (strcmp(str_value, "BOX") == 0 && bounds_verif)
        {
            error = (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);
            for (int i = 0 ; i < 3 ; i++)
            {
                    error = fscanf(conf_file_ptr, "%lf %lf", &(bounds[i][0]), &(bounds[i][1]));
                    if (error != 2)
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
            error = (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);
            for (int i = 0 ; i < *N_atoms ; i++)
            {
                error = fscanf(conf_file_ptr, "%*d %s", str_value);
                if (error != 1)
                {
                    fclose(conf_file_ptr);
                    perror("Atom labels");
                    return EIO;
                }
                error = (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);

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
    int n = 0;
    *pairs = (char **) malloc(sizeof(char **));

    if (*pairs == NULL)
    {
        perror("Pairs first attempt");
        return ENOMEM;
    }

    for (int l1 = 0 ; l1 < N_labels ; l1++)
        for (int l2 = l1 ; l2 < N_labels ; l2++)
        {
            (*pairs)[n] = (char *) malloc(strlen(labels[l1]) + strlen(labels[l2]) + 1);

            if ((*pairs)[n] == NULL)
            {
                perror("Pairs on later attempt");
                return ECANCELED;
            }

            strcpy((*pairs)[n], labels[l1]);
            strcat((*pairs)[n], labels[l2]);
            n++;
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
    double **RDFs = (double **) malloc(N_pairs * sizeof(double *));
    for (int p = 0 ; p < N_pairs ; p++)
    {
        histograms[p] = (int *) calloc(arguments.N_bins, sizeof(int));
        RDFs[p] = (double *) malloc(arguments.N_bins * sizeof(double));
    }

    // The steps-dependent
    double *MSD = (double *) calloc(arguments.N_maxsteps, sizeof(double));
    double *RMSD = (double *) malloc(arguments.N_maxsteps * sizeof(double));

    // The atoms-dependent
    double **initial_upositions = (double **) malloc(N_atoms * sizeof(double *));
    double **current_upositions = (double **) malloc(N_atoms * sizeof(double *));
    double **current_positions = (double **) malloc(N_atoms * sizeof(double *));
    for (int a = 0 ; a < N_atoms ; a++)
    {
        initial_upositions[a] = (double *) malloc(3 * sizeof(double));
        current_upositions[a] = (double *) malloc(3 * sizeof(double));
        current_positions[a] = (double *) malloc(3 * sizeof(double));
    }


    /* Building the hash-tables */
    // Declaring and allocating the arrays
    int **indices = (int **) malloc(N_labels * sizeof(int *));
    int ***pair_indices = (int ***) malloc(N_pairs * sizeof(int **));
    int *pair_numbers = (int *) malloc(N_pairs * sizeof(int));
    bool *pair_status = (bool *) calloc(N_pairs, sizeof(bool));
    int *n_labels = (int *) calloc(N_labels, sizeof(int));
    for (int l = 0 ; l < N_labels ; l++)
        indices[l] = (int *) malloc(N_atoms * sizeof(int));

    // Opening the configuration file
    FILE* conf_file_ptr = fopen(arguments.args[0], "r");
    char dump[STR_BUFF_LIMIT];

    // Error handling
    if (conf_file_ptr == NULL)
    {
        perror("Reading the configuration file to build the hash-tables");
        goto PAIRS;
    }

    // Dumping the first lines
    for (int l = 0 ; l < 9 ; l++)
        err += (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);

    // Reading the first configuration
    for (int a = 0 ; a < N_atoms ; a++)
    {
        // Declaring the temporary variables
        int atom_index;
        char label[STR_LABELS_LIMIT];

        // Reading the first two fields of a line
        err = fscanf(conf_file_ptr, "%d %s", &atom_index, label);
        atom_index--;

        for (int l = 0 ; l < N_labels ; l++)
            if (strcmp(labels[l], label) == 0)
            {
                indices[l][n_labels[l]] = atom_index;
                n_labels[l]++;
            }
        

        /* Reading the initial unwrapped positions */
        err = fscanf(conf_file_ptr, "%*f %*f %*f %lf %lf %lf", &(initial_upositions[atom_index][0]),
                                                                  &(initial_upositions[atom_index][1]),
                                                                  &(initial_upositions[atom_index][2]));
    }

    // Building the pair indices array
    for (int l1 = 0, p = 0; l1 < N_labels ; l1++)
    {
        // indices[l] = (int *) realloc(indices[l], n_labels[l] * sizeof(int));
        for (int l2 = l1 ; l2 < N_labels ; l2++, p++)
        {
            pair_status[p] = (l1 == l2);
            pair_numbers[p] = (!pair_status[p]) ? n_labels[l1] * n_labels[l2] : n_labels[l1] * (n_labels[l1] - 1);
            pair_indices[p] = (int **) malloc(pair_numbers[p] * sizeof(int *));
            
            for (int i1 = 0, j = 0 ; i1 < n_labels[l1] ; i1++)
            {
                for (int i2 = 0 ; i2 < n_labels[l2] ; i2++)
                {
                    if (!pair_status[p] || i1 != i2)
                    {
                        pair_indices[p][j] = (int *) malloc(2 * sizeof(int));
                        pair_indices[p][j][0] = indices[l1][i1],
                        pair_indices[p][j][1] = indices[l2][i2];
                        j++;
                    }
                }
            }
        }
    }

    free(n_labels);
    free(indices);


    /* Processing the configurations */
    // Opening the configuration file
    conf_file_ptr = freopen(arguments.args[0], "r", conf_file_ptr);

    // Error handling
    if (conf_file_ptr == NULL)
    {
        perror("Reading the configuration file to process it");
        goto PAIRS;
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
            err += (fgets(dump, STR_BUFF_LIMIT, conf_file_ptr) == NULL);
        
        // Reading the current configuration
        for (int a = 0 ; a < N_atoms ; a++)
        {
            err = fscanf(conf_file_ptr, "%d %*s", &atom_index);
            atom_index--;
            err = fscanf(conf_file_ptr, "%lf %lf %lf %lf %lf %lf",
                    &(current_positions[atom_index][0]),
                    &(current_positions[atom_index][1]),
                    &(current_positions[atom_index][2]),
                    &(current_upositions[atom_index][0]),
                    &(current_upositions[atom_index][1]),
                    &(current_upositions[atom_index][2]));
            
            // Computing the current MSD
            double diff[3];
            for (int d = 0 ; d < 3 ; d++)
            {
                diff[d] = current_upositions[a][d] - initial_upositions[a][d];
                MSD[step] += diff[d] * diff[d];
            }
        }
        MSD[step] /= (double) N_atoms;

        // Incrementing the histograms
        for (int p = 0 ; p < N_pairs ; p++)
        {
            for (int i = 0 ; i < pair_numbers[p] ; i++)
            {
                double diff[3], distance = 0.;
                for (int d = 0 ; d < 3 ; d++)
                {
                    diff[d] = current_positions[pair_indices[p][i][1]][d] - current_positions[pair_indices[p][i][0]][d];
                    if (diff[d] < -hbounds[d][1])
                        diff[d] += bounds[d][1];
                    else if (diff[d] > hbounds[d][1])
                        diff[d] -= bounds[d][1];
                    distance += diff[d] * diff[d];
                }
                distance = sqrt(distance);
                int index = (int) (distance / delta);
                if (index < arguments.N_bins)
                    histograms[p][index]++;
            }
        }
    }

    // Computing the RMSD for each configuration
    for (int step = 0 ; step < arguments.N_maxsteps ; step++)
    {
        RMSD[step] = (double) sqrt(MSD[step]) / arguments.N_maxsteps;
        MSD[step] /= (double) arguments.N_maxsteps;
    }

    // Freeing the pointers
    free(current_positions);
    free(current_upositions);
    free(initial_upositions);
    fclose(conf_file_ptr);

    // Averaging the RDFs
    // Computing the distance range
    double *r = (double *) malloc(arguments.N_bins * sizeof(double));
    for (int b = 0 ; b < arguments.N_bins ; b++)
        r[b] = (double) b * delta;

    // Normalizing and averaging the histograms
    const double constant = 1. / (4. / 3. * M_PI / V * arguments.N_maxsteps);
    for (int p = 0 ; p < N_pairs ; p++)
    {
        const double pair_constant = constant / pair_numbers[p];

        for (int b = 0 ; b < arguments.N_bins ; b++)
            RDFs[p][b] = pair_constant * histograms[p][b] / (pow(r[b] + delta, 3) - pow(r[b], 3));
    }

    // Shifting the distance range
    for (int b = 0 ; b < arguments.N_bins ; b++)
        r[b] += delta / 2.;


    /* Writing the results */
    FILE *rdfs_file_ptr, *msd_file_ptr;

    // Writing the RDFs
    rdfs_file_ptr = fopen("rdf.dat", "w");

    // Error handling
    if (rdfs_file_ptr == NULL)
    {
        perror("Writing the RDFs");
        goto POINTERS;
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

    // Writing the MSD and RMSD
    msd_file_ptr = fopen("msd.dat", "w");

    // Error handling
    if (msd_file_ptr == NULL)
    {
        perror("Writing the MSD and RMSD");
        goto POINTERS;
    }

    fprintf(msd_file_ptr, "# %7s %7s %7s\n", "step", "MSD", "RMSD");

    for (int step = 0 ; step < arguments.N_maxsteps ; step++)
        fprintf(msd_file_ptr, "  %7d %2.4lf %2.5lf\n", step, MSD[step], RMSD[step]);
    
    fclose(msd_file_ptr);


    /* Exiting normally */
    free(r);
    free(RMSD);
    free(MSD);
    free(RDFs);
    free(histograms);
    free(pairs);
    free(labels);
    exit(EXIT_SUCCESS);


    /* Freeing the pointers and exiting */
    POINTERS:
    free(RMSD);
    free(MSD);
    free(RDFs);
    PAIRS: free(pairs);
    LABELS: free(labels);
    ERROR: errno = err ;
    exit(EXIT_FAILURE);
}