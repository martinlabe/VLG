#include <stdio.h>
#include <igraph/igraph.h>
#include <float.h>
#include <sys/time.h>

#define VERBOSE

igraph_real_t min(igraph_real_t a, igraph_real_t b) {
    igraph_real_t c = a;
    if (a > b)
        c = b;
    return c;
}

igraph_real_t max(igraph_real_t a, igraph_real_t b) {
    igraph_real_t c = a;
    if (a < b)
        c = b;
    return c;
}

igraph_real_t mean(igraph_vector_t *v) {
    igraph_real_t sum = 0;
    for (igraph_integer_t i = 0; i < igraph_vector_size(v); i++)
        sum += VECTOR(*v)[i];
    return sum / igraph_vector_size(v);
}

igraph_real_t var(igraph_vector_t *v, igraph_real_t m) {
    igraph_real_t sum = 0;
    for (igraph_integer_t i = 0; i < igraph_vector_size(v); i++)
        sum += VECTOR(*v)[i] * VECTOR(*v)[i];
    return (sum / igraph_vector_size(v)) - (m * m);
}


igraph_real_t max_vec(igraph_vector_t *v) {
    igraph_real_t m = VECTOR(*v)[0];
    for (igraph_integer_t i = 1; i < igraph_vector_size(v); i++)
        if (m < VECTOR(*v)[i])
            m = VECTOR(*v)[i];
    return m;
}


void read_file(char *filename, igraph_t *graph) {
    FILE *gfile;
    gfile = fopen(filename, "r");
    if (igraph_read_graph_edgelist(graph, gfile, 0, IGRAPH_UNDIRECTED))
    {
#ifdef VERBOSE
        printf("[KO] Reading\n");
#endif
        exit(2);
    }
    fclose(gfile);
#ifdef VERBOSE
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

void get_subset(igraph_vector_t *vertices, igraph_t *graph) {
    /// return a random get_subset of vertices
    igraph_integer_t n, k;
    n = igraph_vcount(graph);
    k = (igraph_integer_t)((double) n * 0.5);

    igraph_vector_init_seq(vertices, 0, igraph_vcount(graph) - 1);
    igraph_vector_shuffle(vertices);
    igraph_vector_resize(vertices, k);
#ifdef VERBOSE
    printf("[OK] Subset of %li elements\n", igraph_vector_size(vertices));
#endif
}

void merge(igraph_t *graph, igraph_vector_t *parents, igraph_vector_t *vdeg) {
    /// XXX uses a lot of RAM by allocationg a vedges
    igraph_vector_t vedges;
    igraph_vector_init(&vedges, igraph_vector_size(parents) * 2);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    long int it = 0;
    for (igraph_integer_t i = 0; i < igraph_vector_size(parents); i++) {
        igraph_integer_t p_i = VECTOR(*parents)[i];
        if (i == p_i)
            continue;
        igraph_neighbors(graph, &neighbors, i, IGRAPH_ALL);

        //get the node degree metric
        VECTOR(*vdeg)[i] = igraph_vector_size(&neighbors);
        unsigned stop = 0;
        for (long int j = 0; j < igraph_vector_size(&neighbors); j++) {
            if (p_i == VECTOR(neighbors)[j]) {
                stop = 1;
                break;
            }
        }
        if (!stop) {
            VECTOR(vedges)[it] = i;
            it += 1;
            VECTOR(vedges)[it] = p_i;
            it += 1;
            //add the new edge in the degree count
            VECTOR(*vdeg)[i]++;
        }
    }
    igraph_vector_resize(&vedges, it);
    igraph_add_edges(graph, &vedges, 0);
    igraph_vector_destroy(&vedges);
}

unsigned nbfs = 0;


int main(int argc, char *argv[]) {

    if (argc != 2)
        exit(1);

    struct timeval t1, t2;
    gettimeofday (&t1, NULL);

    int true_ecc = 8;

    igraph_t graph;
    igraph_t spanner;
    igraph_vector_t subset;

    igraph_vector_t eccmin, eccmax, ecc, eccdelta;
    //vertices degrees
    igraph_vector_t vdeg;

    igraph_vector_t vids, layers, parents;
    igraph_vector_init(&vids, 0);
    igraph_vector_init(&layers, 0);
    igraph_vector_init(&parents, 0);


    read_file(argv[1], &graph);


    igraph_vector_init(&eccmin, igraph_vcount(&graph));
    igraph_vector_fill(&eccmin, -DBL_MAX);
    igraph_vector_init(&eccmax, igraph_vcount(&graph));
    igraph_vector_fill(&eccmax, DBL_MAX);
    igraph_vector_init(&vdeg, igraph_vcount(&graph));
    igraph_vector_init(&ecc, igraph_vcount(&graph));
    igraph_vector_init(&eccdelta, igraph_vcount(&graph));
    igraph_vector_fill(&eccdelta, DBL_MAX);

    // (1) we choose a subset S of get_subset
    get_subset(&subset, &graph);

    // (2) initializing H containing all G get_subset without edges
    igraph_empty(&spanner, igraph_vcount(&graph), IGRAPH_UNDIRECTED);

    // (3) for s in S
    for (igraph_integer_t i = 0; i < igraph_vector_size(&subset); i++) {

        igraph_integer_t v = VECTOR(subset)[i];

        // skip if the ecc have already been found
        if (VECTOR(eccdelta)[i] == 0)
            continue;

        // (3.a) computing the BFS
        igraph_bfs_simple(&graph,
                          v,
                          IGRAPH_ALL,
                          &vids,
                          &layers,
                          &parents);
        nbfs++;

        // (3.b) merging with the spanner
        merge(&spanner, &parents, &vdeg);

        // (3.c) computing metrics between higher and lower bounds
        VECTOR(ecc)[v] = igraph_vector_size(&layers) - 2;
        VECTOR(eccmin)[v] = igraph_vector_size(&layers) - 2;
        VECTOR(eccmax)[v] = igraph_vector_size(&layers) - 2;
        VECTOR(eccdelta)[v] = 0;

        // iterators sur le layer (di = d + 1 without type)
        igraph_integer_t di = 1; // coordinate
        igraph_real_t d = 0; // distance

        // loop on the bfs path
        for (igraph_integer_t j = 0; j < igraph_vector_size(&vids); j++) {

            // the point we are working on
            igraph_integer_t w = VECTOR(vids)[j];

            // compute the distance from the bfs results
            if (j == VECTOR(layers)[di]) {
                di++;
                d++;
            }

            // compute the eccentricity
            if (VECTOR(eccdelta)[w] != 0) {
                VECTOR(eccmin)[w] = max(VECTOR(eccmin)[w], max(VECTOR(ecc)[v] - d, d));
                VECTOR(eccmax)[w] = min(VECTOR(eccmax)[w], VECTOR(ecc)[v] + d);
                VECTOR(eccdelta)[w] = VECTOR(eccmax)[w] - VECTOR(eccmin)[w];
                if (VECTOR(eccdelta)[w] == 0) {
                    VECTOR(ecc)[w] = VECTOR(eccmax)[w];
                }
            }

        }

        if (max_vec(&eccmax) == true_ecc)
        {
            printf("%i\n", i);
            break;
        }
    }

    gettimeofday (&t2, NULL);
    long time = ((t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec) - t1.tv_usec;
    unsigned diameter = (unsigned)max_vec(&eccmax);
    // export the results
    FILE *fp = fopen("../out.txt", "w");
    fprintf(fp, "%u, %u, %ld", diameter, nbfs, time);
    fclose(fp);

    // destroying the bfs objects
    igraph_vector_destroy(&vids);
    igraph_vector_destroy(&layers);
    igraph_vector_destroy(&parents);
    igraph_vector_destroy(&eccmin);
    igraph_vector_destroy(&eccmax);
    igraph_vector_destroy(&ecc);
    igraph_vector_destroy(&eccdelta);
    igraph_vector_destroy(&vdeg);

    // destroy the objects
    igraph_destroy(&graph);
    igraph_destroy(&spanner);

    return 0;
}
