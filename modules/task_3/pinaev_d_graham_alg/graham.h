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

std::vector<point> getRandomArray(size_t size, double max_X = 100.0, 
                                               double max_Y = 100.0);

void Get_Angles(std::vector<double>& ang, std::vector<point>& points);
int Lowest_point(std::vector<point>& points);
void Sort(std::vector<double> &a, std::vector<point>& p);
double result(const point &p1, const point &p2, const point &p3);
bool leftTurn(const point &p1, const point &p2, const point &p3);
void Stack(std::vector<double> &ang, std::vector<point>& points);

//void getConvexHull(std::stack<point> in_stack);
//void getConvexHullParellel(std::stack<point> in_stack);
//bool isConvexHull(std::stack<point> in_stack);

#endif  // MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
