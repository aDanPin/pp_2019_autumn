// Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <math.h>
#include<cmath>
#include<stack>
#include <vector>
#include <random>
#include <algorithm>
#include "../../../modules/task_3/pinaev_d_graham_alg/graham.h"

const double PI = 3.1415;
static int offset = 0;

std::vector<point> getRandomArray(size_t size, int max_X,
                                               int max_Y) {
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);

    // FIXME: resize vector
    std::vector<point> points(size);
    int x, y;
    for(size_t i = 0; i < size; ++i) {
        x = gen() % max_X;
        y = gen() % max_Y;
        points.push_back(point(x, y));
    }

    return points;
}


int LowestPoint(std::vector<point>& points) {
    point first = points[0];
    int first_index = 0;
    int n = points.size();
    for (int i = 1; i < n; ++ i)
        if (points[i].x < first.x ||
            (points[i].x == first.x && points[i].y < first.y)) {
                first = points[i];
                first_index = i;
            }
    return first_index;
}

double area_triangle (point a, point b, point c) {
    return 0.5 * (a.x * b.y + b.x * c.y + c.x * a.y - a.y * b.x - b.y * c.x - c.y * a.x);
}

// ccw - против часововй стрелки
int ccw (point p0, point p1, point p2) {
    return area_triangle(p0, p1, p2) > 0;
}

double dist (point p1, point p2) {
    return std::sqrt(((p1.x - p2.x) * (p1.x - p2.x))
                        + ((p1.x - p2.x) * (p1.x - p2.x)));
}

void Sort(std::vector<point>& p, point first_point) {
    std::sort (p.begin(), p.end()
    , [&](const point& a, const point& b){
        //if (a.index == first.index) return true;
        //if (b.index == first.index) return false;
        if (ccw (first_point, a, b)) return true;
        if (ccw (first_point, b, a)) return false;
        return dist (first_point, a) > dist (first_point, b);
    });
}

std::vector<point> Merge(std::vector<point>& src1, std::vector<point> src2
          , int first, int second
          , point first_point) {

    size_t first_s = first;
    size_t second_s = second;
    if(first_s > src1.size())
        first = src1.size();
    if(second_s > src2.size())
        second = src2.size();

    std::vector<point> dest(first + second);

    int i = 0, j = 0, k = 0;
    while(i < first && j < second) {
        if(ccw (first_point, src1[i], src2[j])){ // src1 is lowest
            dest[k] = src1[i];
            ++i; ++k;
        } else {
            dest[k] = src2[j];
            ++j; ++k;
        }
    }
    // one of src arrs is not completed
    if(i < first) {
        while (i < first) {
            dest[k] = src1[i];
            ++i; ++k;
        }
    } else
    if (j < second) {
        while (j < second) {
            dest[k] = src2[j];
            ++j; ++k;
        }
    }
    return dest;
}

std::vector<point> ParallelSort(std::vector<point>& points, point first_point) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype MPI_Point;
    MPI_Type_contiguous(2, MPI_Point, &MPI_Point);
    MPI_Type_commit(&MPI_Point);

    int n = points.size() / size;
    int res = points.size() % size;
    std::vector<point> dest;

    if (rank == 0) {
        dest.resize(n + res);
    } else {
        dest.resize(n);
    }

    MPI_Scatter(points.data(), n, MPI_Point, dest.data(), n, MPI_Point, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        int a = size - res;
        for (int i = n; i < n + res; ++i) {
            dest[i] = points[a];
            ++a;
        }
        Sort(dest, first_point);
    } else {
        Sort(dest, first_point);
    }

    int locsize = size;
    int i = 0;
    int partner;
    while (locsize > 1) {
        if (locsize % 2 == 0) {
            if (rank % static_cast<int>(pow(2.0, i + 1)) != 0) {
                partner = rank - static_cast<int>(pow(2.0, i));
                MPI_Send(dest.data(), n*static_cast<int>(pow(2.0, i)),
                        MPI_Point, partner, 0, MPI_COMM_WORLD);
                //return std::vector<point>(0);
            } else {
                std::vector<point> tmp(n*static_cast<int>(pow(2.0, i)));
                partner = rank + static_cast<int>(pow(2.0, i));
                MPI_Status status;
                MPI_Recv(tmp.data(), n*static_cast<int>(pow(2.0, i))
                         , MPI_Point, partner, 0, MPI_COMM_WORLD, &status);
                if (rank == 0) {
                    dest = Merge(dest, tmp, dest.size(), n*static_cast<int>(pow(2.0, i))
                                , first_point);
                } else {
                    dest = Merge(dest, tmp, n*static_cast<int>(pow(2.0, i))
                                , n*static_cast<int>(pow(2.0, i))
                                , first_point);
                }
            }
        } else {
            if ((rank % static_cast<int>(pow(2.0, i + 1)) != 0) && (rank != (locsize - 1)*static_cast<int>(pow(2.0, i)))) {
                partner = rank - static_cast<int>(pow(2.0, i));
                MPI_Send(dest.data(), n*static_cast<int>(pow(2.0, i)),
                    MPI_Point, partner, 0, MPI_COMM_WORLD);
                //return std::vector<point>(0);
            } else if (rank == (locsize - 1)*static_cast<int>(pow(2.0, i))) {
                MPI_Send(dest.data(), n*static_cast<int>(pow(2.0, i)),
                    MPI_Point, 0, 0, MPI_COMM_WORLD);
                //return std::vector<point>(0);
            } else {
                if (rank == 0) {
                    std::vector<point> tmp(n*static_cast<int>(pow(2.0, i)));
                    partner = (locsize - 1)*static_cast<int>(pow(2.0, i));
                    MPI_Status status;
                    MPI_Recv(tmp.data(), n*static_cast<int>(pow(2.0, i))
                            , MPI_Point, partner, 0, MPI_COMM_WORLD, &status);
                    dest = Merge(dest, tmp, dest.size(), n*static_cast<int>(pow(2.0, i))
                                , first_point);
                }
                std::vector<point> tmp(n*static_cast<int>(pow(2.0, i)));
                partner = rank + static_cast<int>(pow(2.0, i));
                MPI_Status status;
                MPI_Recv(tmp.data(), n*static_cast<int>(pow(2.0, i))
                        , MPI_Point, partner, 0, MPI_COMM_WORLD, &status);
                if (rank == 0) {
                    dest = Merge(dest, tmp, dest.size(), n*static_cast<int>(pow(2.0, i))
                                , first_point);
                } else {
                    dest = Merge(dest, tmp, n*static_cast<int>(pow(2.0, i))
                                , n*static_cast<int>(pow(2.0, i)), first_point);
                }
            }
        }
        locsize /= 2;
        ++i;
    }

    MPI_Type_free(&MPI_Point);

    return dest;
}

void HullGraham (std::vector<point>& p, std::vector<int> &ip) {
    // первая точка оболочки
    ip.push_back(0);
    // ищем вторую точку оболочки
    double eps = 0.0;
    size_t n = p.size();
    size_t i;
    for(i = 1; i < n && std::abs(area_triangle(p[0], p[1], p[i])) <= eps; ++ i);
    ip.push_back(1);

    // последовательно добавляем точки в оболочку
    int top = 1; // индекс последней точки в оболочке
    while (i < n) {
        // если угол больше pi то извлекаем последнюю точку из оболочки
        if (! ccw (p[ip[top - 1]], p[ip[top]], p[i])) {
            -- top;
            ip.pop_back();
        } else { // иначе кладем ее
            ++ top;
            ip.push_back(i);
            ++ i;
        }
    }

//    for (i = 0; i < ip.size(); ++ i)
//        ip[i] = p[ip[i]].index;
}

bool isConvexHull(std::vector<point>& p, std::vector<int> &ip) {
    if(ip.size() < 3)
        return false;

    bool flag = true;
    for(size_t i = 2; i < ip.size(); ++i)
        if(ccw(p[ip[i]], p[ip[i-1]], p[ip[i-2]])){
            flag = false;
            break;
        }

    return flag;
}

void getConvexHull(std::vector<point>& p, std::vector<int> &ip) {
    int first_index = LowestPoint(p);
    Sort(p, p[first_index]);
    HullGraham(p, ip);
}

void getConvexHullParellel(std::vector<point>& p, std::vector<int> &ip) {
    int first_index = LowestPoint(p);
    p = ParallelSort(p, p[first_index]);
    HullGraham(p, ip);
}
