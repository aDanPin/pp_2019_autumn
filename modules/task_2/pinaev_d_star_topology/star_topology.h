// Copyright 2019 Pinaev Danil
#ifndef MODULES_TASK_2_PINAEV_D_STAR_TOPOLOGY_STAR_TOPOLOGY_H_
#define MODULES_TASK_2_PINAEV_D_STAR_TOPOLOGY_STAR_TOPOLOGY_H_

#include <mpi.h>

MPI_Comm createLinearComm(MPI_Comm old, int size = 0);
bool testLinearTopology(MPI_Comm commLinear);
bool isLinearTopology(MPI_Comm new_comm);

#endif  // MODULES_TASK_2_PINAEV_D_STAR_TOPOLOGY_STAR_TOPOLOGY_H_
