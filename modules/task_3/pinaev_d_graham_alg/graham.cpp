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

const double PI = 3.1415;
static int offset = 0;

void getRandomArray(size_t size, std::vector<point>& vec
                                            , int max_X
                                            , int max_Y) {
    vec.resize(size);
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);

    // FIXME: resize vector
    int x, y;
    for (size_t i = 0; i < vec.size(); ++i) {
        x = gen() % max_X;
        y = gen() % max_Y;
        vec[i] = point(x, y);
    }
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
    // std::cout<<area_triangle(p0, p1, p2)<<"!\n";
    double tmp = area_triangle(p0, p1, p2);
    if (abs(tmp) < 1e-14)
        return 0;
    if (tmp > 0)
        return 1;
    return -1;
    // return area_triangle(p0, p1, p2) >= 0; // >= or > ?????
}

double dist(point p1, point p2) {
    return std::sqrt(((p1.x - p2.x) * (p1.x - p2.x))
                        + ((p1.y - p2.y) * (p1.y - p2.y))); // ЗАМЕНА
}

// void Sort(std::vector<point>& p, size_t first_index) {
//    point first_point = p[first_index];
//    std::swap(p[0], p[first_index]);
//    std::sort (p.begin() + 1, p.end()
//    , [&](const point& a, const point& b){
//        if (ccw (first_point, a, b)) return true;
//        if (ccw (first_point, b, a)) return false;
//        return dist (first_point, a) > dist (first_point, b);
//    });
//    //std::cout<<"Sorted"<<std::endl;
//    //for(size_t i = 0; i < p.size(); ++i)
//    //    std::cout<<p[i].x<<' '<<p[i].y<<std::endl;
// }

void Sort(std::vector<point>& p, point first_point) {
    // check this if std::sort doesn't work:
    // std::vector<point>p1 = p;
    // for(int i = 1; i < p1.size(); i++) {
    //     for(int j = i+1; j < p1.size(); j++) {
    //         if(ccw(first_point, p1[i], p1[j]) < 0)
    //             std::swap(p1[i], p1[j]);
    //     }
    // }
    std::sort(p.begin() + 1, p.end()
    , [&](const point& a, const point& b) {
        return ccw (first_point, a, b) > 0;
        // if (ccw (first_point, a, b)) return true;
        // if (ccw (first_point, b, a)) return false;
        // return dist (first_point, a) <= dist (first_point, b);
    });
    // std::cout<<"Sorted"<<std::endl;
    // for(size_t i = 0; i < p.size(); ++i)
    //    std::cout<<p[i].x<<' '<<p[i].y<<" "<<p1[i].x<<" "<<p1[i].y<<std::endl;
}

std::vector<point> Merge(const std::vector<point>& src1, const std::vector<point>& src2
                        , point first_point) {
    int first = src1.size();
    int second = src2.size();

    std::vector<point> dest(first + second);

    int i = 0, j = 0, k = 0;
    // std::cout<<first<<" "<<second<<"|||\n";
    while (i < first && j < second) {
        // std::cout<<"("<<src1[i].x<<" "<<src1[i].y<<") ("<<src2[j].x<<" "<<src2[j].y<<") "<<i<<", "<<j<<" | "<<ccw(first_point, src1[i], src2[j])<<"\n";
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

void ParallelSort(std::vector<point>& points, point first_point) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype MPI_Point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_Point);
    MPI_Type_commit(&MPI_Point);

/*
    int length, ost;
    if (rank == 0) {
      length = points.size()/size;
      ost = points.size() % size;
    }
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<point> local_vec(length);
    MPI_Scatter(points.data(), length, MPI_Point, local_vec.data(), length, MPI_Point, 0, MPI_COMM_WORLD);
    
    if (rank == 0 && ost != 0) {
        local_vec.insert(local_vec.end(), points.end()-ost, points.end());
    }

    Sort(local_vec, first_point);

    if (rank != 0) {
        MPI_Send(local_vec.data(), length, MPI_Point, 0, 1, MPI_COMM_WORLD);
    } else {
        for (int proc = 1; proc < size; proc++) {
            std::vector<point> tmp_vec(length);
            
            MPI_Status status;
            MPI_Recv(tmp_vec.data(), length, MPI_Point, proc, 1, MPI_COMM_WORLD, &status);
        
            local_vec = Merge(local_vec, tmp_vec, first_point);
        }
    }

    points = local_vec;

    MPI_Type_free (&MPI_Point);
*/
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
                
                break;
                //return local_vec;// ??
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
            
            break;
            //return local_vec;
        }
        iter++;
        i = i/2;
        iter_length *= 2;
    }
    points = local_vec;
}

void HullGraham(std::vector<point>& p, std::vector<int> &ip) {
    // std::cout << "1" << std::endl;
    ip.resize(p.size());
    int* ip_data = ip.data();
    ip_data[0] = 0;
    // std::cout << "2" << std::endl;

    ip_data[1] = 1;

    size_t top = 1;
    // std::cout << "3" << std::endl;
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
    // std::cout << "4" << std::endl;
    ip.resize(top);
}

bool isConvexHull(std::vector<point>& p, std::vector<int> &ip) {
    if (ip.size() < 3)
        return true;

    for (size_t i = 2; i < ip.size(); ++i)
        if (ccw(p[ip[i-2]], p[ip[i-1]], p[ip[i]]) < 0) {
            // std::cout << "Bad points:" << std::endl;
            // std::cout << p[ip[i-2]].x << ' ' << p[ip[i-2]].y << std::endl;
            // std::cout << p[ip[i-1]].x << ' ' << p[ip[i-1]].y << std::endl;
            // std::cout << p[ip[i]].x << ' ' << p[ip[i]].y << std::endl;
            return false;
        }

    return true;
}

bool isSorted(std::vector<point>& p, point first_point) {
    if (p.size() < 3)
        return true;

    for (size_t i = 2; i < p.size(); ++i) {
        if (ccw(first_point, p[i], p[i - 1]) == 1) {
            // std::cout << "On Index " << i << std::endl;

            // std::cout << "Bad points" << std::endl;
            // std::cout << first_point.x << ' ' << first_point.y << std::endl;
            // std::cout << p[i].x << ' ' << p[i].y << std::endl;
            // std::cout << p[i - 1].x << ' ' << p[i - 1].y << std::endl;
            // std::cout << area_triangle(first_point, p[i], p[i-1]) << " " << ccw(first_point, p[i], p[i - 1]) << std::endl;

            return false;
        }
    }
    return true;
}

void getConvexHull(std::vector<point>& p, std::vector<int> &ip) {
    int first_index = LowestPoint(p);
    point first_point = p[first_index];
    point tmp = p[0];
    p[first_index] = tmp;
    p[0] = first_point;
    Sort(p, first_point);
    HullGraham(p, ip);
}
