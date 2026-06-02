#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "data.hpp"
#include <iostream>

/*
enum SolverMethod {
    JACOBI
    SCHWARZ;
}

enum BoundaryConditions {
    HOM_DIRICHLET
    NON_HOM_DIRICHLET
    NEUMANN
    ROBIN
}

enum Mode {
    OPTIMIZED
    MANUAL
}

namespace solver {
template <SolverMethod s, Mode m, BoundaryConditions bc, int Processors = 1, int Threads = 1>
class solver {
    private:
    std::vector<std::vector<double>> create_mesh( ... )
    std::vector<double> set_dirichlet_boundary_conditions( ... )
    std::vector<double> set_neumann_boundary_conditions( ... )
    std::vector<double> set_robin_boundary_conditions( ... )
    std::vector<double> set_boundary_conditions( ... )

    public:
    Solution solve_pde( ... )
    double compute_error( ... )
    void print_solution( ... )
}
*/

Solution solve_pde(ProblemData d, bool par);
void create_vtk(const std::string& filename, const std::vector<double>& u, int n, double h);

#endif