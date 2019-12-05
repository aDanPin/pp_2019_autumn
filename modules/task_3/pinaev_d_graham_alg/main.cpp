//  Copyright 2019 Pinaev Danil
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./graham.h"

TEST(vec_Multiplication, Simple_Multiplication1) {
    vec first(point(0, 0),
                 point(1, 0));

    vec second(point(0, 0),
                  point(0, 1));

    double ans = vecMultiplication(first, second);

    ASSERT_EQ(0.0, ans);
}

TEST(vec_Multiplication, Simple_Multiplication2) {
    vec first(point(1, 2),
                 point(1, 3));

    vec second(point(2, 1),
                  point(3, 1));

    double ans = vecMultiplication(first, second);

    ASSERT_EQ(0.0, ans);
}

TEST(vec_Multiplication, Simple_Multiplication3) {
    vec first(point(2, 1),
                 point(3, 1));

    vec second(point(2, 2),
                  point(3, 2));

    double ans = vecMultiplication(first, second);

    ASSERT_EQ(1.0, ans);
}

TEST(vec_Multiplication, Simple_Multiplication4) {
    vec first(point(2, 1),
                 point(3, 1));

    vec second(point(1, 2),
                  point(1, 1));

    double ans = vecMultiplication(first, second);

    ASSERT_EQ(0.0, ans);
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
