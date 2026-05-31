#include "solver.hpp"
#include<iostream>
#include<vector>
#include<mpi.h>

void run_tests(ProblemData d) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "ne | iter seq | iter par | time seq | time par | speed-up" << std::endl;

    for (int i = 3; i <= 8; ++i) {
        d.ne = std::pow(2, i);
        Solution s1 = solve_pde(d, false);
        Solution s2 = solve_pde(d, true);
        if (rank == 0)
            std::cout << d.ne << " | " << s1.n_iter << " | " << s2.n_iter << " | " << s1.time << " | " << s2.time << " | " << s1.time / s2.time << std::endl; 
    }
}
