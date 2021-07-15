#ifndef UTILS_H
#define UTILS_H

#include <igraph/igraph.h>

igraph_real_t min(igraph_real_t a, igraph_real_t b);
igraph_real_t max(igraph_real_t a, igraph_real_t b);
igraph_real_t mean(igraph_vector_t *v);
igraph_real_t var(igraph_vector_t *v, igraph_real_t m);
igraph_real_t max_vec(igraph_vector_t *v);
igraph_integer_t max_veci(igraph_vector_t *v);
igraph_real_t sum(igraph_vector_t *v);


#endif
