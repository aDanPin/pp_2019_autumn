// Copyright 2019 Pinaev Danil
#ifndef MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_GRAHAM_H_
#define MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_GRAHAM_H_

#include <mpi.h>
#include <vector>

struct point {
    double x, y;

    point(double X, double Y) {
        x = X;
        y = Y;
    }

    point() {
        x = 0;
        y = 0;
    }
};

int LowestPoint(const std::vector<point>& points);  // find lowest point
std::vector<point> Sort(const std::vector<point>& p, point first_point);  // !! // sort points
std::vector<point> ParallelSort(const std::vector<point>& p, point first_point);  // !! // parallel sort points
std::vector<int> HullGraham(const std::vector<point>& p);  // !! // Get a hull
std::vector<point> Merge(const std::vector<point>& src1, const std::vector<point>& src2
                        , point first_point);

int ccw(point p0, point p1, point p2);
double dist(point p1, point p2);
double area_triangle(point a, point b, point c);

std::vector<point> getRandomArray(size_t size, int max_X = 100.0
                                             , int max_Y = 100.0);
bool isConvexHull(const std::vector<point>& p, const std::vector<int> &ip);
bool isSorted(const std::vector<point>& p, point first_point);

#endif  // MODULES_TASK_3_PINAEV_D_GRAHAM_ALG_GRAHAM_H_
