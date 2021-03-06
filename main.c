#include <stdio.h>
#include <stdlib.h>
#include <igraph/igraph.h>
#include <float.h>
#include <sys/time.h>
#include <string.h>
#include "IO.h"
#include "algo.h"
#include "utils.h"




int main(int argc, char *argv[]) {

    if (argc < 3)
        exit(1);

    struct timeval t1, t2;

    igraph_t graph;

    //read the file to the graph
    read_file(argv[2], &graph);


    //get the gcc
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


    //getting the nb of subset in argument
    igraph_integer_t nb_subset = 15;    
    if (argc >= 4)
    {
        nb_subset = atoi(argv[3]);
    }

    // calling the algorithm and taking up the time
    gettimeofday(&t1, NULL);
    struct metrics m = runAlgorithm(atoi(argv[1]), gcc, nb_subset);
    gettimeofday(&t2, NULL);
    long time = ((t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec) - t1.tv_usec;


    // export the results
    char sbuf[1024];
    argv[2][strlen(argv[2]) - 4] = '\0';
    sprintf (sbuf, "%s-%s-%d.txt", argv[1], argv[2] + 10, nb_subset);
    FILE *fp = fopen(sbuf, "a");
    fprintf(fp, "%u, %u, %ld\n", m.diameter, m.nbfs, time);
    fclose(fp);

#if VERBOSE
    printf("\ndiameter = %u\nnb_iter = %u\ntime = %ld\n", m.diameter, m.nbfs, time);
#endif

    // destroy the objects
    igraph_decompose_destroy(&comp);




    return m.diameter;
}
