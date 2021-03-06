#include "utils.h"

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


igraph_real_t sum(igraph_vector_t *v) {
    igraph_real_t s = 0;
    for (igraph_integer_t i = 0; i < igraph_vector_size(v); i++)
        s += VECTOR(*v)[i];
    return s;
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

igraph_integer_t max_veci(igraph_vector_t *v) {
    igraph_integer_t m = 0;
    for (igraph_integer_t i = 1; i < igraph_vector_size(v); i++)
        if (VECTOR(*v)[m] < VECTOR(*v)[i])
            m = i;
    return m;
}

void get_diameter_path(igraph_vector_t *path, igraph_t *spanner, struct bfs_vectors *bfs_v, struct eccentricity *ecc, igraph_integer_t diameter)
{
// can be used to print the path.
    igraph_integer_t max_ecc;
    for (igraph_integer_t j = 0; j < igraph_vector_size(&(ecc->ecc)); j++)
    {
         //get the max eccentricity
         if (VECTOR(ecc->ecc)[j] == diameter)
	 {
	     max_ecc = j;
             break;
	 }
    }
    igraph_bfs_simple(spanner,
                  max_ecc,
                  IGRAPH_ALL,
                  &(bfs_v->vids),
                  &(bfs_v->layers),
                  &(bfs_v->parents));	
    igraph_vector_init(path, 0);
    igraph_integer_t opp = VECTOR(bfs_v->vids)[igraph_vector_size(&bfs_v->vids) - 1];
    igraph_get_shortest_path(spanner, path, NULL, max_ecc,opp, IGRAPH_ALL) ;
    
}



