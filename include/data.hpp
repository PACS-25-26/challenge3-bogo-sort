#ifndef DATA_HPP
#define DATA_HPP

#include <functional>
#include <cmath>
#include <vector>
#include <iostream>

using Vector = std::vector<double>;
using Function = std::function<double(const Vector&)>;


struct ProblemData {

    double x0=0,x1=1,y0=0,y1=1;

    Function f = [](const Vector& x) -> double {
        return 8*M_PI*M_PI*sin(2*M_PI*x[0])*sin(2*M_PI*x[1]);
    };

    Function uex = [](const Vector& x) -> double {
        return sin(2*M_PI*x[0])*sin(2*M_PI*x[1]);
    };

    Function dir_cond = [](const Vector& x) -> double {
        return 0*x[0];
    };

    int ne=128;
    int max_iter=1E6;

    double tol=1E-6;
};


struct Solution {

    Vector u;
    int grid_dimension;
    double h;
    std::vector<std::vector<double>> grid;
    int n_iter;
    double time;

};

#endif