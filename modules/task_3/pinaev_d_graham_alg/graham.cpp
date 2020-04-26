// Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <math.h>
#include <cmath>
#include <stack>
#include <vector>
#include <random>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <cassert>
#include "../../../modules/task_3/pinaev_d_graham_alg/graham.h"

static int offset = 0;

std::vector<point> getRandomArray(size_t size, int max_X, int max_Y) {
    std::vector<point> vec(size);
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);

    int x, y;
    for (size_t i = 0; i < vec.size(); ++i) {
        x = gen() % max_X;
        y = gen() % max_Y;
        vec[i] = point(x, y);
    }
    return vec;
}

int LowestPoint(const std::vector<point>& points) {
    point first = points[0];
    int first_index = 0;
    int n = points.size();
    for (int i = 1; i < n; ++i)
        if (points[i].x < first.x ||
            (points[i].x == first.x && points[i].y < first.y)) {
                first = points[i];
                first_index = i;
            }
    return first_index;
}

double area_triangle(point a, point b, point c) {
    return a.x * b.y + b.x * c.y + c.x * a.y - a.y * b.x - b.y * c.x - c.y * a.x;
}

// ccw - против часовой стрелки
int ccw(point p0, point p1, point p2) {
    double tmp = area_triangle(p0, p1, p2);
    if (std::abs(tmp) < 1e-14)
        return 0;
    if (tmp > 0)
        return 1;
    return -1;
}

double dist(point p1, point p2) {
    return std::sqrt(((p1.x - p2.x) * (p1.x - p2.x))
                        + ((p1.y - p2.y) * (p1.y - p2.y)));
}

std::vector<point> Sort(const std::vector<point>& p, point first_point) {
    std::vector<point> vec(p);
    std::sort(vec.begin(), vec.end()
    , [&](const point& a, const point& b) {
        return ccw (first_point, a, b) > 0;
    });
    return vec;
}

std::vector<point> Merge(const std::vector<point>& src1, const std::vector<point>& src2
                        , point first_point) {
    int first = src1.size();
    int second = src2.size();

    std::vector<point> dest(first + second);

    int i = 0, j = 0, k = 0;
    while (i < first && j < second) {
        if (ccw(first_point, src1[i], src2[j]) >= 0) {  // src1 is lowest
            dest[k] = src1[i];
            ++i; ++k;
        } else {
            dest[k] = src2[j];
            ++j; ++k;
        }
    }
    // one of src arrs is not completed
    while (i < first) {
        dest[k] = src1[i];
        ++i; ++k;
    }
    while (j < second) {
        dest[k] = src2[j];
        ++j; ++k;
    }
    return dest;
}

int pow2(int st) {
  int res = 1;
  for (int i = 0; i < st; i++)
    res *= 2;
  return res;
}

std::vector<point> ParallelSort(const std::vector<point>& points, point first_point) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype MPI_Point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_Point);
    MPI_Type_commit(&MPI_Point);

    int n;
    if (rank == 0) {
        n = points.size() / size;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int>sendcounts(size), displs(size);
    for (int i = 0 ; i < size; i++) {
        sendcounts[i] = points.size()/size + (rank < (signed)(points.size()%size)?1:0);
        if (i != 0)
            displs[i] = displs[i-1] + sendcounts[i-1];
    }
    std::vector<point> dest(sendcounts[rank]);
    MPI_Scatterv(points.data(), sendcounts.data(), displs.data(),
                MPI_Point, dest.data(), sendcounts[rank],
                MPI_Point, 0, MPI_COMM_WORLD);

    dest = Sort(dest, first_point);

    if (rank != 0) {
        MPI_Send(dest.data(), sendcounts[rank], MPI_Point, 0, rank, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        std::vector<point> tmp(n);
        for (int i = 1; i < size; ++i) {
            MPI_Status status;
            MPI_Recv(tmp.data(), sendcounts[i], MPI_Point, MPI_ANY_SOURCE,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            dest = Merge(dest, tmp, first_point);
        }
        dest[0] = first_point;
    }

    MPI_Type_free(&MPI_Point);

    return dest;
}

std::vector<int> HullGraham(const std::vector<point>& p) {
    std::vector<int> ip(p.size());

    int* ip_data = ip.data();
    ip_data[0] = 0;

    ip_data[1] = 1;

    size_t top = 1;
    size_t i = 2;
    size_t n = p.size();
    
    // std::cout<< 1<<std::endl;
    while (i < n) {
        int res = ccw(p[ip[top - 1]], p[ip[top]], p[i]);

//        std::cout<< "Point number (i) "<< i <<std::endl;
//        std::cout<< "Point p[ip[top - 1]] "<< p[ip[top - 1]].x<<' '<<p[ip[top - 1]].y <<std::endl;
//        std::cout<< "Point p[ip[top]] "<<p[ip[top]].x<< ' '<< p[ip[top]].y <<std::endl;
//        std::cout<< "Point p[i] "<< p[i].x << ' '<< p[i].y <<std::endl;
//        std::cout<< "top "<< top <<std::endl;
//        std::cout<< "res "<< res <<std::endl;
//        std::cout <<std::endl <<std::endl;
        

        if (res == 0){ // на одной линии
            ++top;
            ip[top] = i;
            ++i;
        }
        if (res == 1){  // против часовой стрелки
            ++top;
            ip[top] = i;
            ++i;
        }
        if (res == -1){ // по часовой стрелке
            if (top > 1)
                --top;
            else {
                ip[top] = i;
                ++i;
            }

        }


//        if (ccw(p[ip[top - 1]], p[ip[top]], p[i]) < 0) {
//            if (top > 1)
//                --top;
//        } else {
//            ++top;
//            ip[top] = i;
//            ++i;
//        }
    }

//    std::cout<< 2<<std::endl;

    ip.resize(top + 1);
    return ip;
}

bool isConvexHull(const std::vector<point>& p, const std::vector<int> &ip) {
    if (ip.size() < 3)
        return true;

    for (size_t i = 2; i < ip.size(); ++i)
        if (ccw(p[ip[i-2]], p[ip[i-1]], p[ip[i]]) < 0) {
            return false;
        }

    return true;
}

bool isSorted(const std::vector<point>& p, point first_point) {
    if (p.size() < 3)
        return true;

    for (size_t i = 2; i < p.size(); ++i) {
        if (ccw(first_point, p[i], p[i - 1]) == 1) {
            return false;
        }
    }
    return true;
}
