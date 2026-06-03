#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "data.hpp"
#include <iostream>

/*
enum SolverMethod {
    JACOBI
    SCHWARZ;
}


enum Mode {
    OPTIMIZED
    MANUAL
}

namespace solver {
    template <SolverMethod s, Mode m, int Processors = 1, int Threads = 1>
    class pdeSolver {
        private:
            ProblemData d;
        public:

            virtual Solution solve_pde( ... ) = 0;
            double compute_error( ... )
            void create_vtk( ... );
    }
}

*/
double compute_error(ProblemData d, Solution s);
Solution solve_pde(ProblemData d, bool par);
void create_vtk(const std::string& filename, const std::vector<double>& u, int n, double h);

#endif