// Copyright 2019 Pinaev Danil
#ifndef MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
#define MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_

#include <mpi.h>
#include <vector>

struct point{
  double x, y;
  point(double X, double Y){
    x = X;
    y = Y;
  }

  point(){
    x = 0;
    y = 0;
  }
};

std::vector<point> getRandomArray(size_t size, double max_X, double max_Y);
double vectorMultiplication(vec first, vec second);
std::stack<point, std::vector<point>> sortArray(std::vector<point> in_arr);
//void getConvexHull(std::stack<point> in_stack);
//void getConvexHullParellel(std::stack<point> in_stack);
//bool isConvexHull(std::stack<point> in_stack);

#endif  // MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
