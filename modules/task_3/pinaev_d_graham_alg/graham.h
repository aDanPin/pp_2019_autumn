// Copyright 2019 Pinaev Danil
#ifndef MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
#define MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_

#include <mpi.h>
#include <vector>

struct point{
    double x, y;
    //FIXME: index used only for tesing
    point(double X, double Y){
        x = X;
        y = Y;
    }

    point(){
        x = 0;
        y = 0;
    }
};

bool Greater(point a, point b); // point b greater point a
int LowestPoint(std::vector<point>& points); // find lowest point
void Sort(std::vector<point>& p, size_t first_index); // sort points
void Sort(std::vector<point>& p, point first_point); // sort points
void ParallelSort(std::vector<point>& p, size_t first_index); //parallel sort points
void HullGraham (std::vector<point>& p, std::vector<int> &ip); // Get a hull
std::vector<point> Merge(std::vector<point>& src1, std::vector<point> src2
                         , int first, int second
                         , point first_point);

int ccw (point p0, point p1, point p2);
double dist (point p1, point p2);
double area_triangle (point a, point b, point c);

void getRandomArray(size_t size, std::vector<point>& vec
                                , int max_X = 100.0
                                , int max_Y = 100.0);
void getConvexHull(std::vector<point>& p, std::vector<int> &ip);
void getConvexHullParellel(std::vector<point>& p, std::vector<int> &ip);
bool isConvexHull(std::vector<point>& p, std::vector<int> &ip);
bool isSorted(std::vector<point>& p);

#endif  // MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_H_
