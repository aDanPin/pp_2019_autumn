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
    for(size_t i = 0; i < points.size(); ++i) {
        x = gen() % max_X;
        y = gen() % max_Y;
        points[i] = point(x, y);
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
    return area_triangle(p0, p1, p2) >= 0;
}

double dist (point p1, point p2) {
    return std::sqrt(((p1.x - p2.x) * (p1.x - p2.x))
                        + ((p1.x - p2.x) * (p1.x - p2.x)));
}

void Sort(std::vector<point>& p, size_t first_index) {
    point first_point = p[first_index];
    std::swap(p[0], p[first_index]);
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

std::vector<point> ParallelSort(std::vector<point>& points, size_t first_index) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype MPI_Point;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_Point);
    MPI_Type_commit(&MPI_Point);

    int n = points.size() / size;
    int res = points.size() % size;
    std::vector<point> dest;

    dest.resize(n);

    MPI_Scatter(points.data() + res, n, MPI_Point, dest.data(), n, MPI_Point, 0, MPI_COMM_WORLD);

    point first_point = points[first_index];
    Sort(dest, first_index);

    int s = size, n_op = 1;
    while (s > 1)
    {
        s = s/2 + s%2;
        if ((rank-n_op)%(2*n_op) == 0)
        {
            // n - это количество пришедших элементов
            MPI_Send (&n, 1, MPI_INT,
                rank - n_op, 0, MPI_COMM_WORLD);
            // х - это оригинальный массив
            MPI_Send(dest.data(), n, MPI_Point,
                rank - n_op, 0, MPI_COMM_WORLD);
        }
        if ((rank%(2*n_op) == 0) && (size - rank > n_op))
        {
            MPI_Status status;
            int n1;
            MPI_Recv (&n1, 1, MPI_INT,
                rank+n_op, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            std::vector<point> loc_dest(n1);
            MPI_Recv (loc_dest.data(), n1, MPI_Point,
                rank+n_op, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            //for (int i = 0; i < n; i++)
            //    loc_dest[i+n1] = dest[i];
            // bonde - пртосто функция слияния
            dest = Merge(loc_dest, dest
                         , loc_dest.size(), dest.size()
                         , first_point);
            //bond(y, 0, n1-1, n+n1-1);
            //x = new double [n1 + n];
            //for (int i = 0; i < n+n1; i++)
            //    x[i] = y[i];
            n = n + n1;
        }
        n_op*=2;
    }

    MPI_Type_free(&MPI_Point);

    return dest;
}

void HullGraham (std::vector<point>& p, std::vector<int> &ip) {
    std::cout<<"1"<<std::endl;
    // первая точка оболочки
    ip.resize(p.size());
    int* ip_data = ip.data();
    ip_data[0] = 0;
    // ищем вторую точку оболочки
    double eps = 0.0;
    size_t n = p.size();
    size_t i;
    std::cout<<"2"<<std::endl;

    for(i = 1; i < n && std::abs(area_triangle(p[0], p[1], p[i])) <= eps; ++i);
    ip_data[1] = 1;

    // последовательно добавляем точки в оболочку
    size_t top = 1; // индекс последней точки в оболочке
    std::cout<<"3"<<std::endl;
    while (i < n) {
        // если угол больше pi то извлекаем последнюю точку из оболочки
        //std::cout<<"i = "<< i <<std::endl;
        if (! ccw(p[ip[top - 1]], p[ip[top]], p[i])) {
            //std::cout<<"top = "<< top <<" pop. ip.size() = "<<ip.size()<<std::endl;
            -- top;
            //ip.pop_back();
        } else { // иначе кладем ее
            //std::cout<<"top = "<< top <<" push. ip.size() = "<<ip.size()<<std::endl;
            ++ top;
            ip[top] = i;
            ++ i;
        }
    }
    std::cout<<"4"<<std::endl;
    
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

void getConvexHull(std::vector<point>& p, std::vector<int> &ip) {
    int first_index = LowestPoint(p);
    Sort(p, first_index);
    HullGraham(p, ip);
}

void getConvexHullParellel(std::vector<point>& p, std::vector<int> &ip) {
    size_t first_index = LowestPoint(p);
    p = ParallelSort(p, first_index);
    HullGraham(p, ip);
}
