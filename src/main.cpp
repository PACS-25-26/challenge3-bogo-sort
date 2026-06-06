/**
* @file main.cpp
* @brief main file to test the code
*/

#include "data.hpp"
#include "solver.hpp"
#include "tests.hpp"
#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    ProblemData data;
    //run_parallel_vs_serial_tests(data); // uncomment to see the difference between the parallel and the serial version of the code, for different grid dimensions
    run(data);

    MPI_Finalize();
            
}