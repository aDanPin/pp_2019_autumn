//  Copyright 2019 Pinaev Danil
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./graham.h"

const double PI = 3.1415;

TEST(Core_Functionality, Lowest_point) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        // (0, 100) (100, 0) (5, 5) (63, -10) (-10, 0) (0, 0)
        // lowest - is 4
        std::vector<point> points(0);
        points.push_back(point(1000, 100));
        points.push_back(point(100, 0));
        points.push_back(point(5, 5));
        points.push_back(point(63, -100));
        points.push_back(point(-1000, 0));
        points.push_back(point(0, 0));

        int ans = LowestPoint(points);

        ASSERT_EQ(ans, 4);
    }
}

TEST(Core_Functionality, ccw) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::vector<double> ang;
        std::vector<point> points;
        point p22 = point(2, 2);
        point p42 = point(4, 2);
        point p24 = point(2, 4);
        point pm11 = point(-1, 1);
        point p1m1 = point(1, -1);

        ASSERT_EQ(ccw(p22, p24, p42), 0);
        ASSERT_EQ(ccw(p22, p42, p24), 1);
        ASSERT_EQ(ccw(p22, pm11, p1m1), 1);
        ASSERT_EQ(ccw(p22, p1m1, pm11), 0);
    }
}

//TEST(Core_Functionality, Sort) {
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    if (rank == 0) {
//        std::vector<point> points(0);
//
//        point p22 = point(2, 2);
//        point p42 = point(4, 2);
//        point p24 = point(2, 4);
//        point pm11 = point(-1, 1);
//        point p1m1 = point(1, -1);
//
//        points.push_back(p22);
//        points.push_back(p42);
//        points.push_back(p24);
//        points.push_back(pm11);
//        points.push_back(p1m1);
//
//        int first_index = LowestPoint(points);
//        Sort(points, points[first_index]);
//
//        ASSERT_EQ(0, 0); // ??
//        //ASSERT_EQ(points[0].index, pm11.index);
//        //ASSERT_EQ(points[1].index, p1m1.index);
//        //ASSERT_EQ(points[2].index, p42.index);
//        //ASSERT_EQ(points[3].index, p22.index);
//        //ASSERT_EQ(points[4].index, p24.index);
//    }
//}

TEST(Core_Functionality, Merge) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::vector<point> points1(0);
        std::vector<point> points2(0);
        std::vector<point> dest(0);

        points1.push_back(point(0, 0));
        points2.push_back(point(0, 1));
        points1.push_back(point(0, 2));
        points2.push_back(point(0, 3));
        points1.push_back(point(0, 4));
        points2.push_back(point(0, 5));

        point first_point(-1, 0);

        dest = Merge(points1, points2
                    , points1.size(), points2.size()
                    , first_point);

        for(int i = 0; i < 6; ++i)
            ASSERT_EQ(dest[i].y, i);
    }
}

TEST(GrahamAlg, isConvexHull) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::vector<point> points(0);

        point p22 = point(2, 2); // 0
        point p42 = point(4, 2); // 1
        point p24 = point(2, 4); // 2
        point pm11 = point(-1, 1); // 3
        point p1m1 = point(1, -1); // 4

        points.push_back(p22);
        points.push_back(p42);
        points.push_back(p24);
        points.push_back(pm11);
        points.push_back(p1m1);

        // Convex hull - {3, 4, 1, 2}
        std::vector<int>convHull = {3, 4, 1, 2};
        ASSERT_TRUE(isConvexHull(points, convHull));
        
        // Not convex hull - {3, 4, 1, 0, 2}
        std::vector<int>notConvHull = {3, 4, 1, 0, 2};
        ASSERT_FALSE(isConvexHull(points, notConvHull));
    }
}

TEST(Core_Functionality, HullGraham) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::vector<point> points(0);
        std::vector<int> indexes(0);
        point p22 = point(2, 2); // not in hall
        point p42 = point(4, 2); // 3
        point p24 = point(2, 4); // 4
        point pm11 = point(-1, 1); // 1
        point p1m1 = point(1, -1); // 2

        points.push_back(p22);
        points.push_back(p42);
        points.push_back(p24);
        points.push_back(pm11);
        points.push_back(p1m1);

        int first_index = LowestPoint(points);
        Sort(points, first_index);
        HullGraham(points, indexes);

        std::cout<<"Indexes"<<std::endl;
        
        std::cout<<indexes[0]<<std::endl;
        std::cout<<indexes[1]<<std::endl;
        std::cout<<indexes[2]<<std::endl;
        std::cout<<indexes[3]<<std::endl;

        ASSERT_EQ(indexes[0], 0);
        ASSERT_EQ(indexes[1], 1);
        ASSERT_EQ(indexes[2], 2);
        ASSERT_EQ(indexes[3], 4);
    }
}

TEST(GrahamAlg, getConvexHull_Static_Points) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::vector<point> points(0);
        std::vector<int> indexes(0);
        point p22 = point(2, 2); // not in hall
        point p42 = point(4, 2); // 3
        point p24 = point(2, 4); // 4
        point pm11 = point(-1, 1); // 1
        point p1m1 = point(1, -1); // 2

        points.push_back(p22);
        points.push_back(p42);
        points.push_back(p24);
        points.push_back(pm11);
        points.push_back(p1m1);

        getConvexHull(points, indexes);

        ASSERT_TRUE(isConvexHull(points, indexes));
    }
}

TEST(GrahamAlg, getConvexHull_Random_Points) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::vector<point> points = getRandomArray(100);
        std::vector<int> indexes(0);
        std::cout<<"Start Geraham"<<std::endl;
        getConvexHull(points, indexes);
        std::cout<<"Start Check"<<std::endl;
        ASSERT_TRUE(isConvexHull(points, indexes));
    }
}

//TEST(GrahamAlg, getConvexHullParellel_Random_Points) {
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    std::vector<point> points(0);
//    std::vector<int> indexes(0);
//
//    if (rank == 0) {
//        points = getRandomArray(100);
//    }
//
//    getConvexHullParellel(points, indexes);
//
//    if (rank == 0) {
//        ASSERT_TRUE(isConvexHull(points, indexes));
//    }
//}

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
