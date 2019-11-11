// Copyright 2019 Pinaev Danil
#include "../../../modules/task_2/pinaev_d_star_topology/star_topology.h"
#include <mpi.h>

// make star topo comm from old comm
MPI_Comm createStarComm(MPI_Comm old, int size) {
    int nnodes = size;

    int* index = new int[nnodes];
    //first node - N - 1
    for (int i = 0, j = nnodes - 1; i < nnodes; i++, j++)
        index[i] = j;

    int* edges = new int[(nnodes - 1) * 2];
    for (int i = 0; i < nnodes - 1; i++)
        edges[i] = i + 1;
    for (int i = nnodes - 1; i < (nnodes - 1) * 2; i++)
        edges[i] = 0;

    int reorder = 0;
    MPI_Comm commStar;
    MPI_GRAPH_CREATE(old, nnodes, index, edges, reorder, commStar);

    return commStar;
}

bool testStarTopology(MPI_Comm commStar) {
    if (!isStarTopology(commStar))
        return false;
    int rank, size;
    int sourceLess, sourceBig, destLess, destBig;
    int new_coords[LINEARTOPOLEN];
    int val, valFromLess, valFromBig;
    MPI_Status status;
    MPI_Comm_rank(commLinear, &rank);
    MPI_Comm_size(commLinear, &size);
    MPI_Cart_coords(commLinear, rank, LINEARTOPOLEN, new_coords);
    val = new_coords[0];
    valFromLess = -1;
    valFromBig = -1;
    MPI_Cart_shift(commLinear, 0, 1, &sourceLess, &destBig);
    MPI_Cart_shift(commLinear, 0, -1, &sourceBig, &destLess);
    MPI_Sendrecv(&val, 1, MPI_INT, destBig, 4, &valFromLess, 1, MPI_INT,
                 sourceLess, 4, commLinear, &status);
    MPI_Sendrecv(&val, 1, MPI_INT, destLess, 4, &valFromBig, 1, MPI_INT,
                 sourceBig, 4, commLinear, &status);
    MPI_Comm_free(&commLinear);

    if ((rank + 1) != valFromBig) {
        if (new_coords[0] != size - 1) return false;
    }
    if ((rank - 1) != valFromLess) {
        if (new_coords[0] != 0) return false;
    }
    return true;
}

bool isStarTopology(MPI_Comm cur_comm) {
    int topo_type;
    MPI_Topo_test(new_comm, &topo_type);
    if (topology != MPI_GRAPH)
        return false;
    
    int nnodes, nedges;
    MPI_GRAPHDIMS_GET(cur_comm, &nnodes, &nedges);
    
    int index_size = nnodes;
    int edges_size = 2 * (nnodes - 1);
    int* index = new int[index_size];
    int* edges = new int[edges_size];
    MPI_GRAPH_GET(comm, index_size, edges_size, index, edges);
    
    bool flag = true;
    for(int i = 0, j = nnodes - 1; i < index_size; i++, j++)
        if(index[i] != j)
            return false;

    for(int i = 0, j = 1; i < edges_size / 2; i++, j++)
        if(edges[i] != j)
            return false;

    for(int i = nnodes - 1; i < edges_size; i++)
        if(edges[i] != 0)
            return false;

    return true;
}
