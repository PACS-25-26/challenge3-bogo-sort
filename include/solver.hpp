/**
* @file solver.cpp
* @brief declaration file for the PDE solver functions
*/

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "data.hpp"
#include <iostream>

/**
 * @brief This function computes the L2 norm of the error of the solution with respect to the exact solution.
 * @param d data of the problem, containing the domain, the function f, the exact solution uex, the boundary condition function bound_cond, the number of elements ne, the maximum number of iterations max_iter and the tolerance tol.
 * @param s solution of the problem, containing the solution vector u, the grid dimension, the grid itself, the number of iterations, the time taken and the error with respect to the exact solution.
 * @return L2 norm of the error.
 */
double compute_error(ProblemData d, Solution s);

/**
 * @brief This function solves the PDE with dirichlet boundary conditions.
 * @param d data of the problem, containing the domain, the function f, the exact solution uex, the boundary condition function bound_cond, the number of elements ne, the maximum number of iterations max_iter and the tolerance tol.
 * @param par boolean value that indicates whether to further parallelize the code using OpenMP or not.
 * @return solution of the problem, containing the solution vector u, the grid dimension, the grid itself, the number of iterations, the time taken and the error with respect to the exact solution.
*/
Solution solve_pde(ProblemData d, bool par);


/**
 * @brief This function creates a VTK file for visualization.
 * @param filename name of the file to create.
 * @param u solution vector.
 * @param n grid dimension.
 * @param h grid spacing.
 */
void create_vtk(const std::string& filename, const std::vector<double>& u, int n, double h);

#endif