# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>
# include <stdbool.h>


# include "utils.h"


int select_elements(int N_configurations, int *N_atoms, char *labels, Atom **all, int **N_selected, Atom ***selected)
{
    printf("Selecting the elements...\n");


    /* Allocating the array */
    if ((*N_selected = calloc(N_configurations, sizeof(int))) == NULL)
    {
        perror("Allocaitng an array (N_selected)");
        goto NOMEM;
    }

    if ((*selected = malloc(N_configurations * sizeof(Atom *))) == NULL)
    {
        perror("Allocating an array (selected)");
        goto N_SELECTED;
    }

    // Allocating the memory for each configuration
    for (int c = 0 ; c < N_configurations ; c++)
    {
        for (int a = 0 ; a < N_atoms[c] ; a++)
            if (strstr(labels, all[c][a].element) != NULL)
                (*N_selected)[c]++;
        
        if (((*selected)[c] = malloc((*N_selected)[c] * sizeof(Atom))) == NULL)
        {
            perror("Allocating an array slot (selected[])");
            goto SELECTED;
        }
    }


    /* Selecting the atoms */
    for (int c = 0 ; c < N_configurations ; c++)
    {
        printf("conf: %d / %d\r", c + 1, N_configurations);
        int N_selected_local = 0;
        for (int a = 0 ; a < N_atoms[c] && N_selected_local < (*N_selected)[c] ; a++)
        {
            if (strstr(labels, all[c][a].element) == NULL)
                continue;

            if (memcpy(&((*selected)[c][N_selected_local]), &(all[c][a]), sizeof(Atom)) == NULL)
            {
                perror("Copying memory (selected[][], all[][])");
                goto SELECTED_CONF;
            }

            N_selected_local++;
        }
    }

    
    /* Success */
    // Printing the informations
    printf("Selection informations:\n\tSelection: %s\n\tN_selected: %d\n", labels, (*N_selected)[0]);

    // Exiting normally
    return 0;


    /* Errors */
    SELECTED_CONF:
        for (int c = 0 ; c < N_configurations ; c++)
            free(selected[c]);
    SELECTED:
        free(selected);
    N_SELECTED:
        free(N_selected);
    NOMEM:
        return ENOMEM;
}


int select_valency(int N_configurations, int *N_atoms, ComparisonOperator operator, int valency, Atom **all, int **N_selected, Atom ***selected)
{
    printf("Selecting the valency...\n");


    /* Allocating the array */
    if ((*N_selected = calloc(N_configurations, sizeof(int))) == NULL)
    {
        perror("Allocaitng an array (N_selected)");
        goto NOMEM;
    }

    if ((*selected = malloc(N_configurations * sizeof(Atom *))) == NULL)
    {
        perror("Allocating an array (selected)");
        goto N_SELECTED;
    }

    // Allocating the memory for each configuration
    for (int c = 0 ; c < N_configurations ; c++)
    {
        for (int a = 0 ; a < N_atoms[c] ; a++)
        {
            bool condition;
            switch (operator)
            {
                case Greater:
                    condition = (all[c][a].N_bonds > valency);
                    break;
                case GreaterEqual:
                    condition = (all[c][a].N_bonds >= valency);
                    break;
                case Equal:
                    condition = (all[c][a].N_bonds == valency);
                    break;
                case LowerEqual:
                    condition = (all[c][a].N_bonds <= valency);
                    break;
                case Lower:
                    condition = (all[c][a].N_bonds < valency);
                    break;
                default:
                    perror("Selecting the atom");
                    goto N_SELECTED;
            }

            if (condition)
                (*N_selected)[c]++;
        }
        
        if (((*selected)[c] = malloc((*N_selected)[c] * sizeof(Atom))) == NULL)
        {
            perror("Allocating an array slot (selected[])");
            goto N_SELECTED;
        }
    }


    /* Selecting the atoms */
    for (int c = 0 ; c < N_configurations ; c++)
    {
        printf("conf: %d / %d\r", c + 1, N_configurations);
        int N_selected_local = 0;
        for (int a = 0 ; a < N_atoms[c] && N_selected_local < (*N_selected)[c] ; a++)
        {
            bool condition;
            switch (operator)
            {
                case Greater:
                    condition = (all[c][a].N_bonds > valency);
                    break;
                case GreaterEqual:
                    condition = (all[c][a].N_bonds >= valency);
                    break;
                case Equal:
                    condition = (all[c][a].N_bonds == valency);
                    break;
                case LowerEqual:
                    condition = (all[c][a].N_bonds <= valency);
                    break;
                case Lower:
                    condition = (all[c][a].N_bonds < valency);
                    break;
                default:
                    perror("Selecting the atom");
                    goto N_SELECTED;
            }
            if (!condition)
                continue;

            memcpy(&((*selected)[c][N_selected_local]), &(all[c][a]), sizeof(Atom));
            N_selected_local++;
        }
    }

    
    /* Success */
    // Printing the informations
    char *operator_char;
    switch (operator)
    {
        case Greater:
            operator_char = ">";
            break;
        case GreaterEqual:
            operator_char = ">=";
            break;
        case Equal:
            operator_char = "=";
            break;
        case LowerEqual:
            operator_char = "<=";
            break;
        case Lower:
            operator_char = "<";
            break;
    }
    printf("Selection informations:\n\tSelection: valency %s %d\n\tN_selected: %d\n", operator_char, valency, (*N_selected)[0]);

    // Exiting normally
    return 0;


    /* Errors */
    N_SELECTED:
        free(N_selected);
    NOMEM:
        return ENOMEM;
}


int select_coordinate(int N_configurations, int *N_atoms, ComparisonOperator operator, Coordinate coordinate, double value, Atom **all, int **N_selected, Atom ***selected)
{
    printf("Selecting the elements...\n");


    /* Allocating the array */
    if ((*N_selected = calloc(N_configurations, sizeof(int))) == NULL)
    {
        perror("Allocaitng an array (N_selected)");
        return ENOMEM;
    }

    if ((*selected = malloc(N_configurations * sizeof(Atom *))) == NULL)
    {
        perror("Allocating an array (selected)");
        goto N_SELECTED;
    }

    // Allocating the memory for each configuration
    for (int c = 0 ; c < N_configurations ; c++)
    {
        for (int a = 0 ; a < N_atoms[c] ; a++)
        {
            bool condition;
            double atom_coordinate;
            switch (coordinate)
            {
                case X:
                    atom_coordinate = all[c][a].x;
                    break;
                case Y:
                    atom_coordinate = all[c][a].y;
                    break;
                case Z:
                    atom_coordinate = all[c][a].z;
                    break;
                default:
                    perror("Selecting the coordinate");
                    goto N_SELECTED;
            }
            switch (operator)
            {
                case Greater:
                    condition = (atom_coordinate > value);
                    break;
                case GreaterEqual:
                    condition = (atom_coordinate >= value);
                    break;
                case Equal:
                    condition = (atom_coordinate == value);
                    break;
                case LowerEqual:
                    condition = (atom_coordinate <= value);
                    break;
                case Lower:
                    condition = (atom_coordinate < value);
                    break;
                default:
                    perror("Selecting the atom");
                    goto N_SELECTED;
            }

            if (condition)
                (*N_selected)[c]++;
        }
        
        if (((*selected)[c] = malloc((*N_selected)[c] * sizeof(Atom))) == NULL)
        {
            perror("Allocating an array slot (selected[])");
            goto SELECTED;
        }
    }


    /* Selecting the atoms */
    for (int c = 0 ; c < N_configurations ; c++)
    {
        printf("conf: %d / %d\r", c + 1, N_configurations);
        int N_selected_local = 0;
        for (int a = 0 ; a < N_atoms[c] && N_selected_local < (*N_selected)[c] ; a++)
        {
            bool condition;
            double atom_coordinate;
            switch (coordinate)
            {
                case X:
                    atom_coordinate = all[c][a].x;
                    break;
                case Y:
                    atom_coordinate = all[c][a].y;
                    break;
                case Z:
                    atom_coordinate = all[c][a].z;
                    break;
                default:
                    perror("Selecting the coordinate");
                    goto N_SELECTED;
            }

            switch (operator)
            {
                case Greater:
                    condition = (atom_coordinate > value);
                    break;
                case GreaterEqual:
                    condition = (atom_coordinate >= value);
                    break;
                case Equal:
                    condition = (atom_coordinate == value);
                    break;
                case LowerEqual:
                    condition = (atom_coordinate <= value);
                    break;
                case Lower:
                    condition = (atom_coordinate < value);
                    break;
                default:
                    perror("Selecting the atom");
                    goto N_SELECTED;
            }

            if (!condition)
                continue;

            memcpy(&((*selected)[c][N_selected_local]), &(all[c][a]), sizeof(Atom));
            N_selected_local++;
        }
    }

    
    /* Success */
    // Printing the informations
    char *operator_char;
    switch (operator)
    {
        case Greater:
            operator_char = ">";
            break;
        case GreaterEqual:
            operator_char = ">=";
            break;
        case Equal:
            operator_char = "=";
            break;
        case LowerEqual:
            operator_char = "<=";
            break;
        case Lower:
            operator_char = "<";
            break;
    }
    char *coordinate_char;
    switch (coordinate)
    {
        case X:
            coordinate_char = "x";
            break;
        case Y:
            coordinate_char = "y";
            break;
        case Z:
            coordinate_char = "z";
            break;
    }
    printf("Selection informations:\n\tSelection: %s %s %lf\n\tN_selected: %d\n", coordinate_char, operator_char, value, (*N_selected)[0]);

    // Exiting normally
    return 0;


    /* Errors */
    SELECTED:
        free(selected);
    N_SELECTED:
        free(N_selected);
    return ENOMEM;
}


int average_charge(int N_configurations, int *N_atoms, Atom **atoms, Group *group, char *description)
{
    printf("Computing the average charges...\n");


    /* Allocating the arrays */
    if (((*group).N = malloc(N_configurations * sizeof(double))) == NULL)
    {
        perror("Allocating an array (group.N)");
        goto NOMEM;
    }

    if (((*group).average = malloc(N_configurations * sizeof(double))) == NULL)
    {
        perror("Allocating an array (group.average)");
        goto N;
    }


    /* Computing the average charges */
    for (int c = 0 ; c < N_configurations ; c++)
    {
        (*group).average[c] = 0.;
        for (int a = 0 ; a < N_atoms[c] ; a++)
            (*group).average[c] += atoms[c][a].q;
        
        (*group).average[c] /= N_atoms[c];
        (*group).N[c] = N_atoms[c];
    }


    /* Copying the description */
    strncpy((*group).description, description, STR_GROUP_DESCRIPTION_LIMIT - 1);


    /* Success */
    return 0;


    /* Errors */
    N:
        free((*group).N);
    NOMEM:
        return ENOMEM;
}


int write_average(char *file_name, int N_configurations, int *steps, Group group)
{
    printf("Writing the average to '%s'...\n", file_name);


    /* Opening the file */
    FILE *output;
    if ((output = fopen(file_name, "w")) == NULL)
    {
        perror("Opening the file (output)");
        goto IO;
    }


    /* Writing the average */
    fprintf(output, "# %s\n# N_0 = %d\n", group.description, group.N[0]);
    fprintf(output, "# step average\n");
    for (int c = 0 ; c < N_configurations ; c++)
        fprintf(output, "  %d %lf\n", steps[c], group.average[c]);
    

    /* Success */
    fclose(output);

    // Exiting normally
    return 0;


    /* Errors */
    IO:
        return EIO;
}
