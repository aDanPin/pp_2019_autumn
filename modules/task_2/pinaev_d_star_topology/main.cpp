// Copyright 2019 Pinaev Danil
#include <gtest/gtest.h>
#include <gtest-mpi-listener.hpp>
#include <vector>
#include "./star_topology.h"

TEST(Star_Topology_MPI, isStarTopology_Return_Fals_With_Not_Graph_Topo) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_FALSE(isStarTopology(MPI_COMM_WORLD));
    }
}

TEST(Star_Topology_MPI, isStarTopology_Return_Fals_With_Not_Star_Topo) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int* index = new int[size];
    int* edges = new int[size];
    for (int i = 0; i < size; i++) {
        index[i] = i + 1;
        edges[i] = i + 1;
    }

    int reorder = 0;
    MPI_Comm newComm;
    MPI_Graph_create(MPI_COMM_WORLD, size, index, edges, reorder, &newComm);

    if (rank == 0) {
        EXPECT_FALSE(isStarTopology(newComm));
    }
}

TEST(Star_Topology_MPI, isStarTopology_Return_True_With_Star_Topo) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
    MPI_Comm newComm;
    MPI_Graph_create(MPI_COMM_WORLD, nnodes, index, edges, reorder, &newComm);

    if (rank == 0) {
        EXPECT_TRUE(isStarTopology(newComm));
    }
}

TEST(Star_Topology_MPI, Created_Commincator_Is_Star) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm commStar = createStarComm(MPI_COMM_WORLD);
    if (rank == 0) {
        EXPECT_TRUE(isStarTopology(commStar));
    }
}

TEST(Star_Topology_MPI, Created_Communicator_With_Diff_Size_Is_Star) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm commStar = createStarComm(MPI_COMM_WORLD);
    if (rank == 0) {
        EXPECT_TRUE(isStarTopology(commStar));
        MPI_Comm_free(&commStar);
    }
}

TEST(Star_Topology_MPI, Created_Star_Topo_Size_Eq_World_Comm_Size) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm commStar = createStarComm(MPI_COMM_WORLD);
    if (rank == 0) {
        int world_size, star_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_size(commStar, &star_size);

        EXPECT_TRUE(isStarTopology(commStar));
        EXPECT_EQ(world_size, star_size);
        MPI_Comm_free(&commStar);
    }
}

TEST(Star_Topology_MPI, Performance_Test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double stand_topo_start, stand_topo_end, stand_time,
           star_topo_start, star_topo_end, star_time;
    int cur_ans, stand_ans, star_ans;

    char *str;
    size_t strSize = 1000000;
    if (rank == 0) {
        str = getRandomString(strSize);
        cur_ans = getCarNum(str, strSize);
    }


    stand_topo_start = MPI_Wtime();
    stand_ans = getParalCarNum(str, strSize, MPI_COMM_WORLD);
    stand_topo_end = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm commStar = createStarComm(MPI_COMM_WORLD);

    star_topo_start = MPI_Wtime();
    star_ans = getParalCarNum(str, strSize, commStar);
    star_topo_end = MPI_Wtime();

    MPI_Barrier(commStar);
    MPI_Comm_free(&commStar);

    stand_time = stand_topo_end - stand_topo_start;
    star_time = star_topo_end - star_topo_start;

    if (rank == 0) {
        EXPECT_EQ(cur_ans, stand_ans);
        EXPECT_EQ(stand_ans, star_ans);
        EXPECT_GT(stand_time, star_time);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
