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
        points.push_back(point(x, y, i));
    }

    return points;
}


int LowestPoint(std::vector<point>& points) {
    // ищем самую нижнюю из самых левых точек - первая точка
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

void Sort(std::vector<point>& p, int first_index) {
    const point first = p[first_index];
    std::sort (p.begin(), p.end()
    , [&](const point& a, const point& b){
        if (a.index == first.index) return true;
        if (b.index == first.index) return false;
        if (ccw (first, a, b)) return true;
        if (ccw (first, b, a)) return false;
        return dist (first, a) > dist (first, b);
    });
}

void HullGraham (std::vector<point>& p, std::vector<int> &ip) {
    // первая точка оболочки
    ip.push_back(0);
    // ищем вторую точку оболочки
    double eps = 0.0;
    size_t n = p.size();
    size_t i;
    for (i = 1; i < n && std::abs(area_triangle (p[0], p[1], p[i])) <= eps; ++ i);
    ip.push_back (1);

    // последовательно добавляем точки в оболочку
    int top = 1; // индекс последней точки в оболочке
    while (i < n) {
        // если угол больше pi то извлекаем последнюю точку из оболочки
        if (! ccw (p[ip[top - 1]], p[ip[top]], p[i])) {
            -- top;
            ip.pop_back ();
        } else { // иначе кладем ее
            ++ top;
            ip.push_back (i);
            ++ i;
        }
    }
    // ??
    for (i = 0; i < ip.size(); ++ i)
        ip[i] = p[ip[i]].index;
}
