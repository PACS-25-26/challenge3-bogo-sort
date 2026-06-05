/**
* @file tests.hpp
* @brief this file contains the functions used to run tests comparing parallel and serial versions of the PDE solver
* and to run the PDE solver, print results and create visualization files.
*/

#ifndef TESTS_HPP
#define TESTS_HPP

#include "solver.hpp"
#include<iostream>
#include<vector>
#include<mpi.h>

/**
 * @brief This function runs tests comparing parallel and serial versions of the PDE solver.
 * @param d data of the problem, containing the domain, the function f, the exact solution uex, the boundary condition function bound_cond, the number of elements ne, the maximum number of iterations max_iter and the tolerance tol.
 */
void run_parallel_vs_serial_tests(ProblemData d) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "ne |       iter seq |       iter par |       time seq |       time par |       speed-up" << std::endl;

    Solution s1,s2;
    for (int i = 5; i <= 8; ++i) {
        d.ne = std::pow(2, i);
        s1 = solve_pde(d, false);
        s2 = solve_pde(d, true);
        if (rank == 0)
            std::cout << d.ne << " |        " << s1.n_iter << " |       " << s2.n_iter << " |       " << s1.time << " |       " << s2.time << " |       " << s1.time / s2.time << std::endl; 
    }
    
}

/**
 * @brief This function runs the PDE solver, prints results and creates visualization files.
 * @param d data of the problem, containing the domain, the function f, the exact solution uex, the boundary condition function bound_cond, the number of elements ne, the maximum number of iterations max_iter and the tolerance tol.
 */
void run(ProblemData d) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "ne |       iter |       time |" << std::endl;

    Solution s;
    s = solve_pde(d, true);
        if (rank == 0)
            std::cout << d.ne << " |       " << s.n_iter << " |       "<< s.time << " |       " << std::endl; 

    if (rank == 0)
        create_vtk("output/test1.vtk", s.u, s.grid_dimension, s.h);
}

#endif