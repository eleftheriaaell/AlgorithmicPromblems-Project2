#pragma once

#include <vector>
#include <random>
#include <chrono>
#include <bits/stdc++.h>
#include <sys/time.h>
#include "continuous/curve.hpp"
#include "continuous/frechet.hpp"

using namespace std;

struct data{//POINTS NODE
    vector<double> p_data;//VECTOR P
    string ID;
    int d;//DIMENSION
    double current_distance = 0.0;//DISTANCE BETWEEN POINT AND CURRENT QUERY POINT
    int cluster_counter = 0;
    double true_distance = 0.0;
    vector<double> true_curve;
};

bool comp(data const &a, data const &b);

double distance(vector<double> x, vector<double> y, int d);

double Discrete_Frechet_Distance(vector<double> x, vector<double> y);

double ContinuousFrechet_distance(vector<double> current_curve, vector<double> query_curve);

bool compare_distance(data const& a, data const& b);

bool compare_dist(data* const& a, data* const& b);

bool compare_ID_pointer(data* const& a, data* const& b);

//CONTINUOUS 

void filtering(vector<double> current_data, data* current_data_node);