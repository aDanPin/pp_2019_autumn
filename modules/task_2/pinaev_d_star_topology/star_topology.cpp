// Copyright 2019 Pinaev Danil
#include "../../../modules/task_2/pinaev_d_star_topology/star_topology.h"
#include <mpi.h>
#include <random>
#include <string>

// make star topo comm from old comm
MPI_Comm createStarComm(const MPI_Comm old) {
    int size;
    MPI_Comm_size(old, &size);

    int nnodes = size;

    int* index = new int[nnodes];
    for (int i = 0, j = nnodes - 1; i < nnodes; i++, j++)
        index[i] = j;

    int* edges = new int[(nnodes - 1) * 2];
    for (int i = 0; i < nnodes - 1; i++)
        edges[i] = i + 1;
    for (int i = nnodes - 1; i < (nnodes - 1) * 2; i++)
        edges[i] = 0;

    int reorder = 0;
    MPI_Comm commStar;
    MPI_Graph_create(old, nnodes, index, edges, reorder, &commStar);

    return commStar;
}

bool isStarTopology(const MPI_Comm cur_comm) {
    int topo_type;
    MPI_Topo_test(cur_comm, &topo_type);
    if (topo_type != MPI_GRAPH)
        return false;

    int nnodes, nedges;
    MPI_Graphdims_get(cur_comm, &nnodes, &nedges);

    int index_size = nnodes;
    int edges_size = 2 * (nnodes - 1);
    int* index = new int[index_size];
    int* edges = new int[edges_size];
    MPI_Graph_get(cur_comm, index_size, edges_size, index, edges);

    for (int i = 0, j = nnodes - 1; i < index_size; i++, j++)
        if (index[i] != j)
            return false;

    for (int i = 0, j = 1; i < edges_size / 2; i++, j++)
        if (edges[i] != j)
            return false;

    for (int i = nnodes - 1; i < edges_size; i++)
        if (edges[i] != 0)
            return false;

    return true;
}

static int offset = 0;

char* getRandomString(int stringSize) {
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);

    char* str = new char[stringSize];
    std::string charArr = "abcdefghijklmnaoqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890 ";

    for (int i = 0; i < stringSize - 1; ++i)
        str[i] = charArr[gen() % charArr.size()];
    str[stringSize - 1] = '\0';

    return str;
}

int getCarNum(const char* str, int stringSize) {
    int ans = 0;

    for (int i = 0; i < stringSize; ++i)
        if ((str[i] >= 'a' && str[i] <= 'z') ||
            (str[i] >= 'A' && str[i] <= 'Z'))
            ++ans;

    return ans;
}

int getParalCarNum(char* const str, int stringSize, MPI_Comm cur_comm) {
    int size, rank;
    MPI_Comm_size(cur_comm, &size);
    MPI_Comm_rank(cur_comm, &rank);

    int delta;
    int rem;
    if (size > 1) {
        delta = stringSize / (size - 1);
        rem = stringSize % (size - 1);
    } else {
        delta = 0;
        rem = stringSize;
        }
    const char *global_cstr = str;
    char *local_cstr = new char[delta];

    if (rank == 0) {
        for (int proc = 1; proc < size; ++proc)
                MPI_Send(&global_cstr[rem + (proc - 1) * delta], delta, MPI_CHAR, proc, 0, cur_comm);
    } else {
        MPI_Status status;
        MPI_Recv(&local_cstr[0], delta, MPI_CHAR, 0, 0, cur_comm, &status);
    }

    int ans = 0;
    int tmp = 0;
    if (rank == 0) {
        ans = getCarNum(global_cstr, rem);
        for (int proc = 1; proc < size; ++proc) {
            MPI_Status status;
            MPI_Recv(&tmp, 1, MPI_INT, proc, 1, cur_comm, &status);
            ans += tmp;
        }
    } else {
        ans = getCarNum(local_cstr, delta);
        MPI_Send(&ans, 1, MPI_INT, 0, 1, cur_comm);
    }

    delete[] local_cstr;

    return ans;
}
