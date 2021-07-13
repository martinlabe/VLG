#ifndef ALGO_H
#define ALGO_H


#include <igraph/igraph.h>


#define DUMMY 0
#define DELTA_CHOICE 1



struct bfs_vectors
{

    igraph_vector_t vids;
    igraph_vector_t layers;
    igraph_vector_t parents;

};


struct eccentricity
{
    igraph_vector_t eccmin;
    igraph_vector_t eccmax;
    igraph_vector_t ecc;
    igraph_vector_t eccdelta;
};

struct metrics
{

    igraph_integer_t diameter;
    igraph_integer_t nbfs;

};

struct metrics runAlgorithm(int type, igraph_t *graph);

struct metrics algo_dummy(igraph_t *graph, igraph_t *spanner, struct bfs_vectors *bfs_v, struct eccentricity *ecc, igraph_vector_t *vdeg);
struct metrics algo_delta_choice(igraph_t *graph, igraph_t *spanner, struct bfs_vectors *bfs_v, struct eccentricity *ecc, igraph_vector_t *vdeg);


#endif
