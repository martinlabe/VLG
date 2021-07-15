
#include "IO.h"
#include <stdio.h>

void read_file(char *filename, igraph_t *graph) {
    FILE *gfile;
    gfile = fopen(filename, "r");
    if (igraph_read_graph_edgelist(graph, gfile, 0, IGRAPH_UNDIRECTED))
    {
#if VERBOSE
        printf("[KO] Reading\n");
#endif
        exit(2);
    }
    fclose(gfile);
#if VERBOSE
    printf("[OK] Reading %s\n", filename);
    printf("[OK] Graph of %i elements\n", igraph_vcount(graph));
#endif
}

void vector_print(igraph_vector_t *v) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++)
        printf(" %li", (long int) VECTOR(*v)[i]);
    printf("\n");
}

void graph_print(const igraph_t *graph) {
    igraph_vector_t vector;
    igraph_vector_init(&vector, 0);
    igraph_get_edgelist(graph, &vector, 0);
    vector_print(&vector);
}

