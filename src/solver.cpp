/**
* @file solver.cpp
* @brief this file contains the implementation of the functions declared in solver.hpp, 
* which are used to solve the PDE with dirichlet boundary conditions and to create a VTK file for visualization.
*/


#include "solver.hpp"
#include <functional>
#include <cmath>
#include <vector>
#include <iostream>
#include <mpi.h>
#include <fstream>
#include <string>
#include <omp.h>



Solution solve_pde(ProblemData d, bool par) {

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime(); //used to measure the execution time

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Setting the number of threads for OpenMP and making the code serial if par is false.
    if (!par) {
        size = 1;
        omp_set_num_threads(1);
        if (rank != 0) {
            return Solution{};
        }
    }
    else {
        int n_threads = omp_get_num_procs() / size;
        omp_set_num_threads(n_threads);
    }

    int grid_dimension = d.ne;

    double h = (d.x1 - d.x0) / (grid_dimension - 1);

    // global vectors containing the solution at the current and previous iteration, inizialized with 0
    Vector u1(grid_dimension * grid_dimension, 0);
    Vector u0(grid_dimension * grid_dimension, 0);
    
    int iterations = 0;

    double error = d.tol + 1;

    // global grid
    std::vector<std::vector<double>> grid(grid_dimension * grid_dimension);

    // dividing the grid among the processes, handling the reminder 
    // in case the grid dimension is not perfectly divisible by the number of processes
    int proc = grid_dimension / size;
    int remainder = grid_dimension % size;
    std::vector<int> rows_to_rank(size, proc);

    while (remainder > 0) {
        rows_to_rank[remainder] += 1;
        remainder--;
    }
    
    // creating a vector that contains the number of rows assigned to the previous ranks, in order to know from
    //where each rank has to start its computation
    std::vector<int> previous_rows(size, 0);
    for (int i = 1; i < size; ++i) {
        for (int j = i; j < size; ++j) {
            previous_rows[j] += rows_to_rank[i - 1]; 
        }
    }
    
    //initializing variables that will be used in the main loop so that every rank knows which elements it need to skip in the computation
    int local_rows = rows_to_rank[rank];
    int local_previous_rows = previous_rows[rank];
    int skip_start_row = (rank != 0) ? 0 : 1;
    int skip_end_row = (rank != size - 1) ? 0 : 1;

    
    int local_cols = grid_dimension;

    // creating the grid
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < grid_dimension;i++) { 
        for(int j = 0; j < grid_dimension;j++) { 
            grid[i * grid_dimension + j]={d.x0 + i * h, d.y0 + j * h};
        }
    }
    // setting Dirichlet boundary conditions

    if (local_previous_rows == 0) {
        #pragma omp parallel for
        for (int j = 0; j < grid_dimension; ++j)
            u0[j] = d.bound_cond(grid[j]);
        }

    
    if (local_previous_rows + local_rows == grid_dimension) {
        int base = grid_dimension * (grid_dimension - 1);
        #pragma omp parallel for
        for (int j = 0; j < grid_dimension; ++j)
            u0[base + j] = d.bound_cond(grid[base + j]);
    }

    
    #pragma omp parallel for
    for (int i = local_previous_rows; i < local_previous_rows + local_rows; ++i) {
        u0[i * grid_dimension + 0] = d.bound_cond(grid[i * grid_dimension + 0]);
        u0[i * grid_dimension + (grid_dimension - 1)] = d.bound_cond(grid[i * grid_dimension + (grid_dimension - 1)]);
    }



    u1=u0;



    double sum = 0;

    //main loop of the Jacobi method
    #pragma omp parallel
    {
    while(error > d.tol and iterations < d.max_iter) {
        sum = 0;
        

        #pragma omp single 
        {
        //  communicating boundary rows between processes
        if (rank != 0 && par) {
            MPI_Send(&u0[local_previous_rows * grid_dimension], grid_dimension, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u0[local_previous_rows * grid_dimension - grid_dimension], grid_dimension, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        if (rank != size - 1 && par) {
            MPI_Send(&u0[local_previous_rows * grid_dimension + (local_rows - 1) * grid_dimension], grid_dimension, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u0[local_previous_rows * grid_dimension + local_rows * grid_dimension], grid_dimension, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        }

        // updatinh the solution
        #pragma omp for reduction(+:sum) collapse(2) schedule(static)
        for (int i = local_previous_rows + skip_start_row; i < local_previous_rows + local_rows - skip_end_row; ++i) { // for (int i = local_rows * rank + 1; i < local_rows * rank + local_rows - 1; ++i)
            for (int j = 1; j < local_cols - 1; ++j) { // for (j = 1; j < local_cols - 1; ++j)
                u1[i * grid_dimension + j] = (0.25) * (u0[(i + 1) * grid_dimension + j] + u0[(i - 1) * grid_dimension + j] + u0[i * grid_dimension + j + 1] + u0[i * grid_dimension + j - 1] + d.f(grid[i * grid_dimension + j])* h * h);
                sum += std::pow(u1[i * grid_dimension + j] - u0[i * grid_dimension + j],2);
                //e += std::pow(d.f(grid[i * grid_dimension + j]) - u1[i * grid_dimension + j], 2);
            }
        }

        //computing the error and updating the solution for the next iteration
        #pragma omp single 
        {
        error = std::sqrt(h * sum);
        u0 = u1;
        iterations++;

        if (par) {
            MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &iterations, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        }
    }
    }
    
    if (par) {
        MPI_Allreduce(MPI_IN_PLACE, u1.data(), grid_dimension * grid_dimension, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    // measuring the execution time
    if (par)
        MPI_Barrier(MPI_COMM_WORLD);
    
    double t_end = MPI_Wtime();

    double local_time = t_end - t_start;

    if (par)
        MPI_Allreduce(MPI_IN_PLACE, &local_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    
    // assembling the solution in a struct to return it
    Solution s;
    
    s.u = u1;
    s.n_iter = iterations;
    s.grid = grid;
    s.grid_dimension = grid_dimension;
    s.h = h;
    s.time = local_time;
    s.err = compute_error(d, s);

    return s;
}

void create_vtk(const std::string& filename, const std::vector<double>& u, int n, double h) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Errore nell'apertura del file VTK!" << std::endl;
        return;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "Soluzione Equazione di Poisson 2D\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";

    out << "DIMENSIONS " << n << " " << n << " 1\n";
    out << "ORIGIN 0.0 0.0 0.0\n";
    out << "SPACING " << h << " " << h << " 0.0\n";

    out << "POINT_DATA " << n * n << "\n";
    out << "SCALARS u double 1\n";
    out << "LOOKUP_TABLE default\n";

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i * n + j;
            out << u[idx] << "\n";
        }
    }

    out.close();
    std::cout << "File " << filename << " saved!" << std::endl;
}

double compute_error(ProblemData d, Solution s) {
    double error = 0.0;
    for (int i = 0; i < s.grid_dimension; ++i) {
        for (int j = 0; j < s.grid_dimension; ++j) {
            int idx = i * s.grid_dimension + j;
            error += std::pow(d.uex(s.grid[idx]) - s.u[idx], 2);
        }
    }
    return std::sqrt(error);

}