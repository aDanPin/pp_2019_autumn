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

void getRandomArray(size_t size, std::vector<point>& vec
                                            , int max_X
                                            , int max_Y) {
    vec.resize(size);
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);

    // FIXME: resize vector
    int x, y;
    for(size_t i = 0; i < vec.size(); ++i) {
        x = gen() % max_X;
        y = gen() % max_Y;
        vec[i] = point(x, y);
    }
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

//void Sort(std::vector<point>& p, size_t first_index) {
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
//}

void Sort(std::vector<point>& p, point first_point) {
    std::sort (p.begin() + 1, p.end()
    , [&](const point& a, const point& b){
        if (ccw (first_point, a, b)) return true;
        if (ccw (first_point, b, a)) return false;
        return dist (first_point, a) > dist (first_point, b);
    });
    //std::cout<<"Sorted"<<std::endl;
    //for(size_t i = 0; i < p.size(); ++i)
    //    std::cout<<p[i].x<<' '<<p[i].y<<std::endl;  
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
        if(ccw(first_point, src1[i], src2[j])){ // src1 is lowest
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

void ParallelSort(std::vector<point>& points, point first_point) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype MPI_Point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_Point);
    MPI_Type_commit(&MPI_Point);

    int n;
    int res;
    if(rank == 0) {
        n = points.size() / size;
        res = points.size() % size;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //point first_point;
    //if (rank == 0) {
    //    first_point = points[first_index];
    //    point tmp = points[0];
    //    points[0] = first_point;
    //    points[first_index] = tmp;    
    //}
    //MPI_Bcast(&first_point, 1, MPI_Point, 0, MPI_COMM_WORLD);

    std::vector<point> dest;
    if(rank == 0){
        dest.resize(n + res);
    } else {
        dest.resize(n);
    }
   
    MPI_Scatter(points.data() + res, n, MPI_Point, dest.data(), n, MPI_Point, 0, MPI_COMM_WORLD);

    if(rank == 0){
        for(int i = 0; i < n + res; ++i)
            dest[i] = points[i];
    }

    Sort(dest, first_point);

    if(rank != 0){
        MPI_Send(dest.data(), n, MPI_Point, 0, rank, MPI_COMM_WORLD);
    }
    if(rank == 0){
        std::vector<point> tmp(n);
        for(int i = 1; i < size; ++i){
            MPI_Status status;
            MPI_Recv(tmp.data(), n, MPI_Point, MPI_ANY_SOURCE,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            dest = Merge(dest, tmp, dest.size(), tmp.size(), first_point);
        }
        dest[0] = first_point;
    }
    points = dest;
    
    //return points;
}

void HullGraham (std::vector<point>& p, std::vector<int> &ip) {
    std::cout<<"1"<<std::endl;
    ip.resize(p.size());
    int* ip_data = ip.data();
    ip_data[0] = 0;
    std::cout<<"2"<<std::endl;

    ip_data[1] = 1;

    size_t top = 1;
    std::cout<<"3"<<std::endl;
    size_t i = 2;
    size_t n = p.size();
    while (i < n) {
        if (! ccw(p[ip[top - 1]], p[ip[top]], p[i])) {
            -- top;
        } else {
            ++ top;
            ip[top] = i;
            ++ i;
        }
    }
    std::cout<<"4"<<std::endl;
    ip.resize(top);
}

bool isConvexHull(std::vector<point>& p, std::vector<int> &ip) {
    if(ip.size() < 3)
        return true;

    for(size_t i = 2; i < ip.size(); ++i)
        if(!ccw(p[ip[i-2]], p[ip[i-1]], p[ip[i]])){
            std::cout<<"Bad points:"<<std::endl;
            std::cout<<p[ip[i-2]].x<<' '<<p[ip[i-2]].y<<std::endl;
            std::cout<<p[ip[i-1]].x<<' '<<p[ip[i-1]].y<<std::endl;
            std::cout<<p[ip[i]].x<<' '<<p[ip[i]].y<<std::endl;
            return false;
        }

    return true;
}

bool isSorted(std::vector<point>& p, point first_point) {
    if(p.size() < 3)
        return true;

    for(size_t i = 2; i < p.size(); ++i)
        if(ccw(first_point, p[i], p[i - 1])) {
            std::cout<<"On Index "<<i<<std::endl;
            
            std::cout<<"Bad points"<<std::endl;
            std::cout<<first_point.x<<' '<<first_point.y<<std::endl;
            std::cout<<p[i].x<<' '<<p[i].y<<std::endl;
            std::cout<<p[i - 1].x<<' '<<p[i - 1].y<<std::endl;
            std::cout<<area_triangle(first_point, p[i], p[i-1])<<std::endl;
            
            return false;
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
