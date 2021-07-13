#ifndef IO_H
#define IO_H

#include <igraph/igraph.h>

void read_file(char *filename, igraph_t *graph);

void vector_print(igraph_vector_t *v);

void graph_print(const igraph_t *graph);

#endif
