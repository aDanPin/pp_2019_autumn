// Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <math.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_3/pinaev_d_graham_alg/graham.h"

static int offset = 0;

std::vector<point> getRandomArray(size_t size, int max_X, int max_Y){
    std::mt19937 gen;
    gen.seed((unsigned)time(0) + ++offset);

    std::vector<point> points(size);

    for(size_t i = 0; i < size; ++i){
        double x = gen() % max_X;
        double y = gen() % max_Y;
        points[i] = point(x, y);
    }

    return points;
}

double vectorMultiplication(vec first, vec second){
    double first_x = first.second_p.x - first.first_p.x;
    double first_y = first.second_p.y - first.first_p.y;
    double second_x = second.second_p.x - second.first_p.x;
    double second_y = second.second_p.y - second.first_p.y;
    
    double res = first_x * second_x + first_y * second_y;

    return res;
}

std::stack<point, std::vector<point>> sortArray(std::vector<point> in_arr) {
    //find left point
    for(size_t i = 1; i < in_arr.size(); ++i){
        if(in_arr[i].x < in_arr[0].x){
            point tmp = in_arr[0];
            in_arr[0] = in_arr[i];
            in_arr[i] = tmp;
        }
    }

    //sort points
    vec v0(in_arr[0], point(0, 0));
    vec v1, v2;
    //sort vstavkami
    for (size_t i = 1; i < in_arr.size(); ++i) {
        //find max
        point p_max = in_arr[i];
        for (size_t j = i; j < in_arr.size(); ++j) {
            v1 = vec(in_arr[0], p_max);
            v2 = vec(in_arr[0], in_arr[j]);

            if (vectorMultiplication(v0, v2) >
                  vectorMultiplication(v0, v1)) {
                      point tmp = p_max;
                      p_max = in_arr[j];
                      in_arr[j] = tmp;
            }
      }
    }
}
//void getConvexHull(std::stack<point> in_stack);
//void getConvexHullParellel(std::stack<point> in_stack);
//bool isConvexHull(std::stack<point> in_stack);
