#include "data.hpp"
#include "solver.hpp"
#include "tests.hpp"
#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    ProblemData data;
    run_tests(data);

    MPI_Finalize();
            
}