// Copyright 2019 Pinaev Danil
#ifndef MODULES_TASK_2_PINAEV_D_STAR_TOPOLOGY_STAR_TOPOLOGY_H_
#define MODULES_TASK_2_PINAEV_D_STAR_TOPOLOGY_STAR_TOPOLOGY_H_

#include <mpi.h>

MPI_Comm createStarComm(const MPI_Comm old);
bool isStarTopology(const MPI_Comm new_comm);

char* getRandomString(int stringSize);
int getCarNum(const char* str, int stringSize);
int getParalCarNum(char* const str, int stringSize, const MPI_Comm cur_comm);

#endif  // MODULES_TASK_2_PINAEV_D_STAR_TOPOLOGY_STAR_TOPOLOGY_H_
