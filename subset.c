#include <igraph/igraph.h>
#include "IO.h"
#include "utils.h"
#include "subset.h"

//return a completly random subset
void get_subset(igraph_vector_t *vertices, const igraph_t *graph) {
    igraph_integer_t n, k;
    n = igraph_vcount(graph);

    igraph_vector_init_seq(vertices, 0, igraph_vcount(graph) - 1);
    igraph_vector_shuffle(vertices);
}

//return a subset containing a random point per community, with the communites in random order
void get_community_subset(igraph_vector_t *res, const igraph_t *graph) {
    
    //calculate communities
    igraph_vector_t membership;
    igraph_vector_init(&membership, 0);
    igraph_integer_t nb_clusters = 1;
    igraph_real_t gsize = igraph_vcount(graph);
    igraph_real_t div = (0.5/gsize);
    igraph_community_leiden(graph, NULL, NULL, div,0.01, 0, &membership, &nb_clusters, NULL);



    igraph_vector_t vertices;
    igraph_vector_bool_t found;


    igraph_vector_init_seq(&vertices, 0, igraph_vcount(graph) - 1); 
    igraph_vector_init(res, nb_clusters);
    igraph_vector_bool_init(&found, nb_clusters);

    
    igraph_vector_shuffle(&vertices);
    igraph_real_t count = 0;
   

    for (igraph_integer_t i = 0; i < igraph_vcount(graph); i++)
    {
	igraph_integer_t pos = VECTOR(vertices)[i];
	igraph_integer_t c = VECTOR(membership)[pos];
        if(VECTOR(found)[c] == 0)
	{
	    VECTOR(found)[c] = 1;
	    VECTOR(*res)[c] = pos;
	    count++;
	    if (count == nb_clusters)
	        break;
	}
    } 
    igraph_vector_shuffle(res);

    igraph_vector_destroy(&membership);
    igraph_vector_bool_destroy(&found);
}

//return a random point per community subset with communities in descending order relative to their sizes
void get_max_community_subset(igraph_vector_t *res, const igraph_t *graph) {
    

    igraph_vector_t membership;
    igraph_vector_init(&membership, 0);
    igraph_integer_t nb_clusters = 1;
    igraph_real_t gsize = igraph_vcount(graph);
    igraph_real_t div = (0.5/gsize);
    igraph_community_leiden(graph, NULL, NULL, div,0.01, 0, &membership, &nb_clusters, NULL);

    igraph_vector_t vertices, histo, loc;
    igraph_vector_bool_t found;



    igraph_vector_init_seq(&vertices, 0, igraph_vcount(graph) - 1); 
    igraph_vector_init(&histo, nb_clusters);
    igraph_vector_init(&loc, nb_clusters);
    igraph_vector_init(res, nb_clusters);
    igraph_vector_shuffle(&vertices);

    for (igraph_integer_t i = 0; i < igraph_vcount(graph); i++)
    {
	igraph_integer_t pos = VECTOR(vertices)[i];
	igraph_integer_t c = VECTOR(membership)[pos];
	VECTOR(histo)[c]++;
	if (VECTOR(histo)[c] == 1)
	    VECTOR(loc)[c] = pos;
    }
    //fonction de sort peu propre, mais negligeable
    for (igraph_integer_t i = 0; i < igraph_vector_size(&histo); i++)
    {
        igraph_integer_t loc_max = max_veci(&histo);
	VECTOR(*res)[i] = VECTOR(loc)[loc_max];
	VECTOR(histo)[loc_max] = 0;
    }

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&vertices);
    igraph_vector_destroy(&histo);
    igraph_vector_destroy(&loc);
}
