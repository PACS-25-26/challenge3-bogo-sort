/**
* @file scalability_test.cpp
* @brief main file to test the scalability of the parallel code
*/

#include "data.hpp"
#include "solver.hpp"
#include <iostream>
#include <mpi.h>
#include <iomanip>

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ProblemData data;

    Solution s = solve_pde(data, true);

    //std::cout << std::setw(10) << std::setw(15) << std::setw(15) << std::setw(15) << std::setw(15);
    if(rank==0)
     std::cout << "TIME=" << s.time
              << " ERROR=" << s.err
              << " ITER=" << s.n_iter
              << std::endl;

    MPI_Finalize();
            
}