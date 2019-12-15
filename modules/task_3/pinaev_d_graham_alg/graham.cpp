// Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <math.h>
#include <cmath>
#include <stack>
#include <vector>
#include <random>
#include <algorithm>
#include <ctime>
#include <cassert>
#include "../../../modules/task_3/pinaev_d_graham_alg/graham.h"

const double PI = 3.1415;
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
    if (abs(tmp) < 1e-14)
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

    int length, ost;
    if (rank == 0) {
      length = points.size()/size;
      ost = points.size() % size;
    }

    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<point>local_vec(length);

    if (rank == 0 && ost != 0) {
        local_vec.insert(local_vec.end(), points.end()-ost, points.end());
    }

    int i = size;
    int iter = 1;
    int iter_length = length;
    int partner;
    int loc_size = length;
    if (rank == 0)
    loc_size += ost;
    while (i > 1) {
        if (i % 2 == 1) {
            if (rank == 0) {
                std::vector<point>neig_vec(iter_length);

                partner = pow2(iter - 1) * (i - 1);

                MPI_Status status;
                MPI_Recv(neig_vec.data(), iter_length, MPI_Point, partner, 0, MPI_COMM_WORLD, &status);

                local_vec = Merge(local_vec, neig_vec, first_point);

                loc_size += iter_length;
            }
            if (rank == pow2(iter - 1) * (i - 1)) {
                MPI_Send(local_vec.data(), iter_length, MPI_Point, 0, 0, MPI_COMM_WORLD);

                return local_vec;
            }
        }

        if (rank % pow2(iter) == 0) {
            partner = rank + pow2(iter-1);

            std::vector<point> neig_vec(iter_length);

            MPI_Status status;
            MPI_Recv(neig_vec.data(), iter_length, MPI_Point, partner, 1, MPI_COMM_WORLD, &status);

            local_vec = Merge(local_vec, neig_vec, first_point);

            loc_size += iter_length;
        }
        if (rank % pow2(iter) != 0) {
            partner = rank - pow2(iter-1);

            MPI_Send(local_vec.data(), iter_length, MPI_Point, partner, 1, MPI_COMM_WORLD);

            return local_vec;
        }
        iter++;
        i = i/2;
        iter_length *= 2;
    }
    return local_vec;
}

std::vector<int> HullGraham(const std::vector<point>& p) {
    std::vector<int> ip(p.size());

    int* ip_data = ip.data();
    ip_data[0] = 0;

    ip_data[1] = 1;

    size_t top = 1;
    size_t i = 2;
    size_t n = p.size();
    while (i < n) {
        if (ccw(p[ip[top - 1]], p[ip[top]], p[i]) < 0) {
            --top;
        } else {
            ++top;
            ip[top] = i;
            ++i;
        }
    }
    ip.resize(top);
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
