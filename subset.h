#ifndef SUBSET_H
#define SUBSET_H


#include <igraph/igraph.h>

void get_subset(igraph_vector_t *vertices, const igraph_t *graph);

void get_community_subset(igraph_vector_t *res, const igraph_t *graph);

void get_max_community_subset(igraph_vector_t *res, const igraph_t *graph);

#endif
