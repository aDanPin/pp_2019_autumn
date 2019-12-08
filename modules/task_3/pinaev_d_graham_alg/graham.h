// Copyright 2019 Pinaev Danil
#ifndef MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
#define MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_

#include <mpi.h>
#include <vector>

struct point{
    double x, y;
    //FIXME: index used only for tesing
    int index;
    point(double X, double Y, int idx = 0){
        x = X;
        y = Y;
        index = idx;
    }

    point(){
        x = 0;
        y = 0;
        index = 0;
    }
};

bool Greater(point a, point b); // point b greater point a
int LowestPoint(std::vector<point>& points); // find lowest point
void Sort(std::vector<point>& p, int first_index); // sort points
void ParallelSort(std::vector<point>& p, int first_index); //parallel sort points
void HullGraham (std::vector<point>& p, std::vector<int> &ip); // Get a hull

int ccw (point p0, point p1, point p2);
double dist (point p1, point p2);
double area_triangle (point a, point b, point c);

std::vector<point> getRandomArray(size_t size, int max_X = 100.0,
                                               int max_Y = 100.0);
void getConvexHull(std::vector<point>& p, std::vector<int> &ip);
//void getConvexHullParellel(std::stack<point> in_stack);
bool isConvexHull(std::vector<point>& p, std::vector<int> &ip);

#endif  // MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
