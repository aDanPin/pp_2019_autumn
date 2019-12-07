//  Copyright 2019 Pinaev Danil
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./graham.h"

const double PI = 3.1415;
//
//TEST(Core_Functionality, Lowest_point) {
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    if (rank == 0) {
//        // (0, 100) (100, 0) (5, 5) (63, -10) (-10, 0) (0, 0)
//        // lowest - is 4
//        std::vector<point> points;
//        points.push_back(point(0, 100));
//        points.push_back(point(100, 0));
//        points.push_back(point(5, 5));
//        points.push_back(point(63, -100));
//        points.push_back(point(-10, 0));
//        points.push_back(point(0, 0));
//
//        int ans = Lowest_point(points);
//
//        ASSERT_EQ(ans, 3);
//    }
//}

//TEST(Core_Functionality, Get_Angles) {
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    // FIXME: check accuracy  
//    if (rank == 0) {
//        std::vector<double> ang;
//        std::vector<point> points;
//        points.push_back(point(-1, 0));
//        points.push_back(point(3, -4));
//        points.push_back(point(2, -3));
//        points.push_back(point(1, -2));
//
//        //int ans = Lowest_point(points);
//        //ASSERT_EQ(ans, 0);
//
//        Get_Angles(ang, points);
//        ASSERT_NEAR(135.003, ang[0], 0.001);
//        ASSERT_NEAR(0.0, ang[1], 0.001);
//        ASSERT_NEAR(135.003, ang[2], 0.001);
//        ASSERT_NEAR(135.003, ang[3], 0.001);
//    }
//}

//TEST(Core_Functionality, Sort) {
//    ASSERT_EQ(0.0, ans);
//}
//
//TEST(Core_Functionality, Result) {
//    ASSERT_EQ(0.0, ans);
//}
//
//TEST(Core_Functionality, leftTurn) {
//    ASSERT_EQ(0.0, ans);
//}
//
//TEST(Core_Functionality, Stack) {
//    ASSERT_EQ(0.0, ans);
//}
//
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
