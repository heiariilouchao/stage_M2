# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>


# include "utils.h"


int select_elements(int N_configurations, int *N_atoms, char *labels, Atom **all, int **N_selected, Atom ***selected)
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
            if (strstr(labels, all[c][a].element) != NULL)
                (*N_selected)[c]++;
        
        if (((*selected)[c] = malloc((*N_selected)[c] * sizeof(Atom))) == NULL)
        {
            perror("Allocating an array slot (selected[])");
            goto NOMEM;
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
                goto NOMEM;
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
    N_SELECTED: free(N_selected);
    NOMEM: free(selected);
    return ENOMEM;
}
