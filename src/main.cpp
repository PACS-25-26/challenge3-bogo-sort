#include "data.hpp"
#include "solver.hpp"
#include <iostream>

int main(int argc, char* argv[]){
    ProblemData data;
    Solution s = solve_pde(data, argc, argv);

       // if(s.n_iter<data.max_iter) std::cout << "The method has converged after "<< s.n_iter << " iterations" << std::endl;

    create_vtk("output/test1.vtk", s.u, s.grid_dimension, s.h);
            
}