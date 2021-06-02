#include <stdio.h>
#include <igraph/igraph.h>

void vector_print(igraph_vector_t *v)
{
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}

void graph_print(const igraph_t *graph)
{
    igraph_vector_t vector;
    igraph_vector_init(&vector, 0);
    igraph_get_edgelist(graph, &vector, 0);
    vector_print(&vector);
}

void merge(igraph_t *graph, igraph_vector_t *parents)
{
    /// XXX uses a lot of RAM by allocationg a vedges
    igraph_vector_t vedges;
    igraph_vector_init(&vedges, igraph_vector_size(parents) * 2);

    long int it = 0;
    for (long int i = 0; i < igraph_vector_size(parents); i++)
    {
        long int p_i = VECTOR(*parents)[i];
        igraph_vector_t neighbors;
        igraph_vector_init(&neighbors, 0);
        igraph_neighbors(graph, &neighbors, i, IGRAPH_ALL);
        unsigned stop = 0;
        if (i == p_i)
            continue;
        for (long int j = 0; j < igraph_vector_size(&neighbors); j++)
        {
            if (p_i == VECTOR(neighbors)[j])
            {
                stop = 1;
                break;
            }
        }
        if (!stop)
        {
            VECTOR(vedges)[it] = i;
            it += 1;
            VECTOR(vedges)[it] = p_i;
            it += 1;
        }
    }
    igraph_vector_resize(&vedges, it);
    igraph_add_edges(graph, &vedges, 0);
    igraph_vector_destroy(&vedges);
}


int main() {
    printf("Hello, World!\n");

    igraph_t graph;
    igraph_t spanner;

    // reading the file
    FILE *gfile;
    gfile = fopen("/home/martin/Data/VLG/graph_easy.txt","r");
    if(igraph_read_graph_edgelist(&graph, gfile, 0, IGRAPH_UNDIRECTED))
        printf("rip\n");

    // initializing variables
    igraph_empty(&spanner, igraph_vcount(&graph), IGRAPH_UNDIRECTED);

    // bfs
    igraph_vector_t vids, layers, parents;
    igraph_vector_init(&vids, 0);
    igraph_vector_init(&layers, 0);
    igraph_vector_init(&parents, 0);

    igraph_bfs_simple(&graph,
                      5,
                      IGRAPH_ALL,
                      &vids,
                      &layers,
                      &parents);

    printf("BFS:\n");
    vector_print(&vids);
    vector_print(&layers);
    vector_print(&parents);
    printf("\n");

    // merge
    printf("MERGE:");
    merge(&spanner, &parents);
    graph_print(&spanner);

    // destroy the objects
    igraph_destroy(&graph);
    igraph_destroy(&spanner);
    igraph_vector_destroy(&vids);
    igraph_vector_destroy(&layers);
    igraph_vector_destroy(&parents);

    fclose(gfile);
    return 0;
}
