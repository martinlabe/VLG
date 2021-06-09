#include <stdio.h>
#include <igraph/igraph.h>

void read_file(igraph_t *graph)
{
    FILE *gfile;
    gfile = fopen("/home/martin/Data/VLG/graph_easy.txt","r");
    if(igraph_read_graph_edgelist(graph, gfile, 0, IGRAPH_UNDIRECTED))
        printf("rip\n");
    fclose(gfile);
}

void vector_print(igraph_vector_t *v)
{
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++)
        printf(" %li", (long int) VECTOR(*v)[i]);
    printf("\n");
}

void graph_print(const igraph_t *graph)
{
    igraph_vector_t vector;
    igraph_vector_init(&vector, 0);
    igraph_get_edgelist(graph, &vector, 0);
    vector_print(&vector);
}

void get_subset(igraph_vector_t *vertices, igraph_t *graph)
{
    /// return a random get_subset of vertices
    igraph_integer_t n, k;
    n = igraph_vcount(graph);
    k = (igraph_integer_t)((double)n * 0.5);

    igraph_vector_init_seq(vertices, 0, igraph_vcount(graph) - 1);
    igraph_vector_shuffle(vertices);
    igraph_vector_resize(vertices, k);
}

void merge(igraph_t *graph, igraph_vector_t *parents)
{
    /// XXX uses a lot of RAM by allocationg a vedges
    igraph_vector_t vedges;
    igraph_vector_init(&vedges, igraph_vector_size(parents) * 2);

    long int it = 0;
    for (igraph_integer_t i = 0; i < igraph_vector_size(parents); i++)
    {
        igraph_integer_t p_i = VECTOR(*parents)[i];
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
    igraph_vector_t subset;

    igraph_vector_t vids, layers, parents;
    igraph_vector_init(&vids, 0);
    igraph_vector_init(&layers, 0);
    igraph_vector_init(&parents, 0);

    read_file(&graph);

    // (1) we choose a subset S of get_subset
    printf("SUBSET:\n");
    get_subset(&subset, &graph);
    vector_print(&subset);

    // (2) initializing H containing all G get_subset without edges
    igraph_empty(&spanner, igraph_vcount(&graph), IGRAPH_UNDIRECTED);

    // (3) for s in S
    for (igraph_integer_t i = 0; i < igraph_vector_size(&subset); i++)
    {
        // (3.a) computing the BFS
        igraph_bfs_simple(&graph,
                          VECTOR(subset)[i],
                          IGRAPH_ALL,
                          &vids,
                          &layers,
                          &parents);

        // (3.b) merging with the spanner
        merge(&spanner, &parents);

        // (3.c) computing the difference between higher and lower bounds
    }

    graph_print(&spanner);

    // destroying the bfs objects
    igraph_vector_destroy(&vids);
    igraph_vector_destroy(&layers);
    igraph_vector_destroy(&parents);
    // destroy the objects
    igraph_destroy(&graph);
    igraph_destroy(&spanner);

    return 0;
}
