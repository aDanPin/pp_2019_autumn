//  Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./graham.h"

#include <iostream>

const double PI = 3.1415;

TEST(GrahamAlg, getConvexHullParellel_Random_Points) {
    int rank;
    //MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//    int size = atoi(argv[1]);

    std::vector<point> points;
    std::vector<int> indexes;

    size_t first_index;
    point first_point;
    if (rank == 0) {
        points = getRandomArray(100);

        first_index = LowestPoint(points);
        point first_point = points[first_index];
        point tmp = points[0];
        points[first_index] = tmp;
        points[0] = first_point;
    }

    points = ParallelSort(points, first_point); // parralel shit

    if (rank == 0) {
        // Sort checking
        bool isSorted = true;
        if (points.size() < 3)
            isSorted = true;
        else {
            for (size_t i = 2; i < points.size(); ++i) {
                if (ccw(first_point, points[i], points[i - 1]) == 1) {
                    isSorted = false;
                    break;
                }
            }
        }

        if(isSorted)
            std::cout<< "Array is sorted"<<std::endl;
        else 
            std::cout<< "Array is not sorted!!!"<<std::endl;


//      HULL GRAHAM
        indexes.resize(points.size());

        int* ip_data = indexes.data();
        ip_data[0] = 0;

        ip_data[1] = 1;

        size_t top = 1;
        size_t i = 2;
        size_t n = points.size();
    
        while (i < n) {
            int res = ccw(points[indexes[top - 1]], points[indexes[top]], points[i]);

            if (res == 0){ // на одной линии
                ++top;
                indexes[top] = i;
                ++i;
            }
            if (res == 1){  // против часовой стрелки
                ++top;
                indexes[top] = i;
                ++i;
            }
            if (res == -1){ // по часовой стрелке
                if (top > 1)
                    --top;
                else {
                    indexes[top] = i;
                    ++i;
                }

            }
        }

        indexes.resize(top + 1);

        bool isConvex = true;
        if (indexes.size() < 3)
            isConvex = true;
        else {
            for (size_t i = 2; i < indexes.size(); ++i)
                if (ccw(points[indexes[i-2]], points[indexes[i-1]], points[indexes[i]]) < 0) {
                    isConvex = false;
                    break;
                }
        }

        if (isConvex)
            std::cout<< "Convex is hull" << std::endl;
        else 
            std::cout<< "Convex is not hull!!!" << std::endl;
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
