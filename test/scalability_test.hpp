#include "data.hpp"
#include "solver.hpp"
#include <iostream>
#include <mpi.h>
#include <iomanip>

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    ProblemData data;
    data.ne = 100;

    Solution s = solve_pde(data, true);

    std::cout << std::setw(10) << std::setw(15) << std::setw(15) << std::setw(15) << std::setw(15);

    MPI_Finalize();
            
}