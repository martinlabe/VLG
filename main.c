#include <stdio.h>
#include <igraph/igraph.h>

void print_vector(const igraph_vector_t* vids)
{
    unsigned n = igraph_vector_size(vids);
    for (unsigned i = 0; i < n; i++)
    {
        printf("%u ", (unsigned)VECTOR(*vids)[i]);
    }
    printf("\n");
}

int main() {
    printf("Hello, World!\n");

    // initializing variables
    igraph_t graph;

    FILE *gfile;
    gfile = fopen("graph_easy.txt","r");



    if(!igraph_read_graph_edgelist(graph, gfile, 0, IGRAPH_UNDIRECTED))
        printf("rip\n")

    igraph_vector_t vids;
    igraph_vector_init(&vids, 0);
    igraph_vector_t parents;
    igraph_vector_init(&parents, 0);

    igraph_bfs_simple(&graph,
                      1,
                      0,
                      &vids,
                      0,
                      &parents);

    print_vector(&vids);
    print_vector(&parents);

    // destroy the objects

    fclose(gfile)
    igraph_destroy(&graph);
    igraph_vector_destroy(&vids);
    return 0;
}
