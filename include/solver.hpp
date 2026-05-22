#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "data.hpp"
#include <iostream>

Solution solve_pde(ProblemData d);
void create_vtk(const std::string& filename, const std::vector<double>& u, int n, double h);

#endif