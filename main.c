#include <stdio.h>
#include <stdlib.h>
#include <igraph/igraph.h>
#include <float.h>
#include <sys/time.h>
#include "IO.h"
#include "algo.h"
#define VERBOSE




int main(int argc, char *argv[]) {

    if (argc != 3)
        exit(1);

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);


    igraph_t graph;

    read_file(argv[2], &graph);


    struct metrics m = runAlgorithm(atoi(argv[1]), &graph);

    gettimeofday(&t2, NULL);
    long time = ((t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec) - t1.tv_usec;
    // export the results
    FILE *fp = fopen("../out.txt", "w");
    fprintf(fp, "%u, %u, %ld", m.diameter, m.nbfs, time);
    fclose(fp);

    // destroy the objects
    igraph_destroy(&graph);

    return 0;
}
