#include <stdio.h>
#include <stdlib.h>
#include <igraph/igraph.h>
#include <float.h>
#include <sys/time.h>
#include "IO.h"
#include "algo.h"
#include "utils.h"
#define VERBOSE




int main(int argc, char *argv[]) {

    if (argc < 3)
        exit(1);

    struct timeval t1, t2;


    igraph_t graph;



    //read the file to the graph
    read_file(argv[2], &graph);







    igraph_vector_ptr_t comp;
    igraph_vector_ptr_init(&comp, 0);
    
    igraph_decompose(&graph, &comp, IGRAPH_STRONG, -1, 3);

    igraph_destroy(&graph);


    igraph_t *gcc = VECTOR(comp)[0];

    for (igraph_integer_t i = 1; i < igraph_vector_ptr_size(&comp); i++)
    {
        if (igraph_vcount(gcc) < igraph_vcount(VECTOR(comp)[i]))
	    gcc = VECTOR(comp)[i];
    }


/*    igraph_vector_t test;
    igraph_vector_init(&test, 0);
    igraph_vs_t testv;
    igraph_vs_1(&testv, 1125203);
    igraph_eccentricity(gcc, &test, testv, IGRAPH_ALL);

    printf("test0 = %f", VECTOR(test)[0]);

    return 0;*/

    igraph_integer_t nb_subset = 15;

    if (argc >= 4)
    {
        nb_subset = atoi(argv[3]);
    }
    
    gettimeofday(&t1, NULL);
    
     
    struct metrics m = runAlgorithm(atoi(argv[1]), gcc, nb_subset);
    
    
    gettimeofday(&t2, NULL);
    long time = ((t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec) - t1.tv_usec;


    // export the results
    char sbuf[1024];
    sprintf (sbuf, "out_%s", argv[2] + 10);
    FILE *fp = fopen(sbuf, "w");
    fprintf(fp, "%u, %u, %ld\n", m.diameter, m.nbfs, time);
    fclose(fp);

    printf("%u, %u, %ld\n", m.diameter, m.nbfs, time);
   



    // destroy the objects
    igraph_decompose_destroy(&comp);




    return 0;
}
