#include "data.hpp"
#include "solver.hpp"
#include <iostream>

int main(){
    ProblemData data;
    Solution s= solve_pde(data);

    if(s.n_iter<data.max_iter) std::cout << "The method has converged after "<< s.n_iter << " iterations" << std::endl;
    create_vtk("output/test1.vtk", s.u, s.grid_dimension, s.h);
            
}