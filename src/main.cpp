#include "data.hpp"
#include "solver.hpp"
#include "tests.hpp"
#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    ProblemData data;
    //Solution s = solve_pde(data, false);

    // if(s.n_iter<data.max_iter) std::cout << "The method has converged after "<< s.n_iter << " iterations" << std::endl;

    //create_vtk("output/test1.vtk", s.u, s.grid_dimension, s.h);
    run_tests(data);

    MPI_Finalize();
            
}