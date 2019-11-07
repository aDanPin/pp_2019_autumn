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

bool testLinearTopology(MPI_Comm commLinear) {
    if (!isLinearTopology(commLinear)) return false;
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

bool isLinearTopology(MPI_Comm new_comm) {
    int curDims, topology;
    MPI_Topo_test(new_comm, &topology);
    if (topology != MPI_CART) return false;
    MPI_Cartdim_get(new_comm, &curDims);
    if (curDims != LINEARTOPOLEN) return false;

    int dims[LINEARTOPOLEN], periods[LINEARTOPOLEN], coords[LINEARTOPOLEN];
    MPI_Cart_get(new_comm, LINEARTOPOLEN, dims, periods, coords);
    if (periods[0] != 0) return false;
    return true;
}
