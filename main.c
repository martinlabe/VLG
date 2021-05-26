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
    igraph_vector_t v;
    igraph_real_t edges[] = {1, 2, 1, 3, 1, 4, 1, 6, 2, 5, 2, 3, 2, 6, 3, 4, 4, 5, 4, 6};

    // create a vector carrying more efficiently the edges
    igraph_vector_view(&v,
                       edges,
                       sizeof(edges) / sizeof(double));

    // create the graph
    igraph_create(&graph,
                  &v,
                  0,
                  IGRAPH_UNDIRECTED);

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
    igraph_destroy(&graph);
    igraph_vector_destroy(&vids);
    return 0;
}
