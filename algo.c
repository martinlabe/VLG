#include <igraph/igraph.h>
#include <float.h>
#include "algo.h"
#include "IO.h"
#include "utils.h"

void get_subset(igraph_vector_t *vertices, const igraph_t *graph) {
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


struct metrics runAlgorithm(int type, igraph_t *graph) {


//    int true_ecc = 8;

    igraph_t spanner;


    igraph_vector_t vids, layers, parents;

    //vectors used by bfs
    struct bfs_vectors bfs_v;

    igraph_vector_init(&(bfs_v.vids), 0);
    igraph_vector_init(&(bfs_v.layers), 0);
    igraph_vector_init(&(bfs_v.parents), 0);


    struct eccentricity ecc;
    
    
    //vertices degrees
    igraph_vector_t vdeg;

    igraph_vector_init(&(ecc.eccmin), igraph_vcount(graph));
    igraph_vector_fill(&(ecc.eccmin), -DBL_MAX);

    igraph_vector_init(&(ecc.eccmax), igraph_vcount(graph));
    igraph_vector_fill(&(ecc.eccmax), DBL_MAX);

    igraph_vector_init(&vdeg, igraph_vcount(graph));

    igraph_vector_init(&(ecc.ecc), igraph_vcount(graph));

    igraph_vector_init(&(ecc.eccdelta), igraph_vcount(graph));
    igraph_vector_fill(&(ecc.eccdelta), DBL_MAX);


    //run selected algo
    
    struct metrics m;

    if (type == DUMMY)
        m = algo_dummy(graph, &spanner, &bfs_v, &ecc, &vdeg);
    else if (type == DELTA_CHOICE)
        m = algo_delta_choice(graph, &spanner, &bfs_v, &ecc, &vdeg);
 
    // destroying the bfs objects
    igraph_vector_destroy(&(bfs_v.vids));
    igraph_vector_destroy(&(bfs_v.layers));
    igraph_vector_destroy(&(bfs_v.parents));
    igraph_vector_destroy(&(ecc.eccmin));
    igraph_vector_destroy(&(ecc.eccmax));
    igraph_vector_destroy(&(ecc.ecc));
    igraph_vector_destroy(&(ecc.eccdelta));
    igraph_vector_destroy(&vdeg);

    // destroy the objects
    igraph_destroy(&spanner);

    return m;
}

struct metrics algo_dummy(igraph_t *graph, igraph_t *spanner, struct bfs_vectors *bfs_v, struct eccentricity *ecc, igraph_vector_t *vdeg)
{
    struct metrics m = {0,0};

    // (1) we choose a random subset S of get_subset
    igraph_vector_t subset;
    
    get_subset(&subset, graph);

    // (2) initializing H containing all G get_subset without edges
    igraph_empty(spanner, igraph_vcount(graph), IGRAPH_UNDIRECTED);

    // (3) for s in S
    for (igraph_integer_t i = 0; i < igraph_vector_size(&subset); i++) {

        igraph_integer_t v = VECTOR(subset)[i];

        // skip if the ecc have already been found
        if (VECTOR(ecc->eccdelta)[i] == 0)
            continue;

        // (3.a) computing the BFS
        igraph_bfs_simple(graph,
                          v,
                          IGRAPH_ALL,
                          &(bfs_v->vids),
                          &(bfs_v->layers),
                          &(bfs_v->parents));
        m.nbfs++;

        // (3.b) merging with the spanner
        merge(spanner, &(bfs_v->parents), vdeg);

        // (3.c) computing metrics between higher and lower bounds
        VECTOR(ecc->ecc)[v] = igraph_vector_size(&(bfs_v->layers)) - 2;
        VECTOR(ecc->eccmin)[v] = igraph_vector_size(&(bfs_v->layers)) - 2;
        VECTOR(ecc->eccmax)[v] = igraph_vector_size(&(bfs_v->layers)) - 2;
        VECTOR(ecc->eccdelta)[v] = 0;

        // iterators sur le layer
        igraph_integer_t d = 0; // distance

        // loop on the bfs path
        for (igraph_integer_t j = 0; j < igraph_vector_size(&(bfs_v->vids)); j++) {

            // the point we are working on
            igraph_integer_t w = VECTOR(bfs_v->vids)[j];

            // compute the distance from the bfs results
            if (j == VECTOR(bfs_v->layers)[d + 1]) {
                d++;
            }

            // compute the eccentricity
            if (VECTOR(ecc->eccdelta)[w] != 0) {
                VECTOR(ecc->eccmin)[w] = max(VECTOR(ecc->eccmin)[w], max(VECTOR(ecc->ecc)[v] - d, d));
                VECTOR(ecc->eccmax)[w] = min(VECTOR(ecc->eccmax)[w], VECTOR(ecc->ecc)[v] + d);
                VECTOR(ecc->eccdelta)[w] = VECTOR(ecc->eccmax)[w] - VECTOR(ecc->eccmin)[w];
                if (VECTOR(ecc->eccdelta)[w] == 0) {
                    VECTOR(ecc->ecc)[w] = VECTOR(ecc->eccmax)[w];
                }
            }

        }
        m.diameter = max_vec(&(ecc->eccmax)); 
        if (m.diameter == 8) //to replace with a good function
        {
            printf("%i\n", i);
            break;
        }
    }

    igraph_vector_destroy(&subset);
    return m;

}

struct metrics algo_delta_choice(igraph_t *graph, igraph_t *spanner, struct bfs_vectors *bfs_v, struct eccentricity *ecc, igraph_vector_t *vdeg)
{
    struct metrics m = {0,0};

    // (1) we choose a random subset S of get_subset
    
    //faire un vrai random à la place
    igraph_vector_t subset;
    get_subset(&subset, graph);
    igraph_integer_t next = VECTOR(subset)[0];



    // (2) initializing H containing all G get_subset without edges
    igraph_empty(spanner, igraph_vcount(graph), IGRAPH_UNDIRECTED);

    // (3) for s in S //mettre une vrai condition d'arret
    for (igraph_integer_t i = 0; i < igraph_vector_size(&subset); i++) {

        igraph_integer_t v = next;

        // skip if the ecc have already been found
        if (VECTOR(ecc->eccdelta)[i] == 0)
            continue;

        // (3.a) computing the BFS
        igraph_bfs_simple(graph,
                          v,
                          IGRAPH_ALL,
                          &(bfs_v->vids),
                          &(bfs_v->layers),
                          &(bfs_v->parents));
        m.nbfs++;

        // (3.b) merging with the spanner
        merge(spanner, &(bfs_v->parents), vdeg);

        // (3.c) computing metrics between higher and lower bounds
        VECTOR(ecc->ecc)[v] = igraph_vector_size(&(bfs_v->layers)) - 2;
        VECTOR(ecc->eccmin)[v] = igraph_vector_size(&(bfs_v->layers)) - 2;
        VECTOR(ecc->eccmax)[v] = igraph_vector_size(&(bfs_v->layers)) - 2;
        VECTOR(ecc->eccdelta)[v] = 0;

        // iterators sur le layer
        igraph_integer_t d = 0; // distance

        // loop on the bfs path
        for (igraph_integer_t j = 0; j < igraph_vector_size(&(bfs_v->vids)); j++) {

            // the point we are working on
            igraph_integer_t w = VECTOR(bfs_v->vids)[j];

            // compute the distance from the bfs results
            if (j == VECTOR(bfs_v->layers)[d + 1]) {
                d++;
            }

            // compute the eccentricity
            if (VECTOR(ecc->eccdelta)[w] != 0) {
                VECTOR(ecc->eccmin)[w] = max(VECTOR(ecc->eccmin)[w], max(VECTOR(ecc->ecc)[v] - d, d));
                VECTOR(ecc->eccmax)[w] = min(VECTOR(ecc->eccmax)[w], VECTOR(ecc->ecc)[v] + d);
                VECTOR(ecc->eccdelta)[w] = VECTOR(ecc->eccmax)[w] - VECTOR(ecc->eccmin)[w];

                if (VECTOR(ecc->eccdelta)[w] > VECTOR(ecc->eccdelta)[next])
		    next = w;



		if (VECTOR(ecc->eccdelta)[w] == 0) {
                    VECTOR(ecc->ecc)[w] = VECTOR(ecc->eccmax)[w];
                }
            }

        }
        m.diameter = max_vec(&(ecc->eccmax)); 
        if (m.diameter == 8) //to replace with a good function
        {
            printf("%i\n", i);
            break;
        }
    }

    igraph_vector_destroy(&subset);
    return m;

}


