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
    double t_start = MPI_Wtime();

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (!par) {
        size = 1;
        omp_set_num_threads(1);
        if (rank != 0) {
            return Solution{};
        }
    }
    else { // && Mode == OPTIMIZED
        int n_threads = omp_get_num_procs() / size;
        omp_set_num_threads(n_threads);
    }

    int grid_dimension = d.ne;

    double h = (d.x1 - d.x0) / (grid_dimension - 1);

    // questi saranno i vettori global
    Vector u1(grid_dimension * grid_dimension, 0);
    Vector u0(grid_dimension * grid_dimension, 0);
    
    int iterations = 0;

    // questo sarà l'errore globale
    double error = d.tol + 1;

    // questa sarà la griglia globale
    std::vector<std::vector<double>> grid(grid_dimension * grid_dimension);

    // per implementare la divisione delle righe tra i processi anche se le righe (e quindi le colonne) sono dispari, si usa la strategia:
    int proc = grid_dimension / size;
    int remainder = grid_dimension % size;
    std::vector<int> rows_to_rank(size, proc);

    while (remainder > 0) {
        rows_to_rank[remainder] += 1;
        remainder--;
    }
    
    std::vector<int> previous_rows(size, 0);
    for (int i = 1; i < size; ++i) {
        for (int j = i; j < size; ++j) {
            previous_rows[j] += rows_to_rank[i - 1]; 
        }
    }
    
    int local_rows = rows_to_rank[rank];
    int local_previous_rows = previous_rows[rank];
    int skip_start_row = (rank != 0) ? 0 : 1;
    int skip_end_row = (rank != size - 1) ? 0 : 1;

    // in questo modo creo: 
    // - rows_to_rank, ovvero un vettore contenente il numero di righe destinate ad ogni rank;
    // - previous_rows, ovvero un vettore che contiene in posizione i il numero di righe precedenti al rank i 
    // - skip_start_row e skip_end_row, ovvero due valori che servono al primo rank per farlo partire una riga dopo e all'ultimo per farlo fermare una riga prima 

    int local_cols = grid_dimension;

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < grid_dimension;i++) { // sostituire grid_dimension con local_rows
        for(int j = 0; j < grid_dimension;j++) { // sostituire grid_dimension con local_cols (== grid_dimension)
            grid[i * grid_dimension + j]={d.x0 + i * h, d.y0 + j * h};
        }
    }

    if (local_previous_rows == 0) {
        #pragma omp parallel for
        for (int j = 0; j < grid_dimension; ++j)
            u0[j] = d.bound_cond(grid[j]);
        }

    // bottom row (ultima riga globale) se il rank possiede l'ultima riga
    if (local_previous_rows + local_rows == grid_dimension) {
        int base = grid_dimension * (grid_dimension - 1);
        #pragma omp parallel for
        for (int j = 0; j < grid_dimension; ++j)
            u0[base + j] = d.bound_cond(grid[base + j]);
    }

    // prime/ultime colonne per le righe locali
    #pragma omp parallel for
    for (int i = local_previous_rows; i < local_previous_rows + local_rows; ++i) {
        u0[i * grid_dimension + 0] = d.bound_cond(grid[i * grid_dimension + 0]);
        u0[i * grid_dimension + (grid_dimension - 1)] = d.bound_cond(grid[i * grid_dimension + (grid_dimension - 1)]);
    }

    u1=u0;



    double sum = 0;
    //double e = 0;

    #pragma omp parallel
    {
    while(error > d.tol and iterations < d.max_iter) {
        sum = 0;
        //e = 0;

        #pragma omp single 
        {
        //  il rank 0 non ha rank precedenti
        if (rank != 0 && par) {
        //  devo inviare la prima riga di u1 del mio rank al rank precedente
            MPI_Send(&u0[local_previous_rows * grid_dimension], grid_dimension, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        //  devo ricevere l'ultima riga di u1 del rank precedente e memorizzarla nella riga di u1 precedente alla prima riga occupata
            MPI_Recv(&u0[local_previous_rows * grid_dimension - grid_dimension], grid_dimension, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // il rank massimo non ha rank successivi
        if (rank != size - 1 && par) {
        //  devo inviare l'ultima riga di u1 del mio rank al rank successivo
            MPI_Send(&u0[local_previous_rows * grid_dimension + (local_rows - 1) * grid_dimension], grid_dimension, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        //  devo ricevere la prima riga di u1 del rank successivo e memorizzarla nella riga di u1 successiva all'ultima riga occupata
            MPI_Recv(&u0[local_previous_rows * grid_dimension + local_rows * grid_dimension], grid_dimension, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        }

        #pragma omp for reduction(+:sum) collapse(2) schedule(static)
        for (int i = local_previous_rows + skip_start_row; i < local_previous_rows + local_rows - skip_end_row; ++i) { // for (int i = local_rows * rank + 1; i < local_rows * rank + local_rows - 1; ++i)
            for (int j = 1; j < local_cols - 1; ++j) { // for (j = 1; j < local_cols - 1; ++j)
                u1[i * grid_dimension + j] = (0.25) * (u0[(i + 1) * grid_dimension + j] + u0[(i - 1) * grid_dimension + j] + u0[i * grid_dimension + j + 1] + u0[i * grid_dimension + j - 1] + d.f(grid[i * grid_dimension + j])* h * h);
                sum += std::pow(u1[i * grid_dimension + j] - u0[i * grid_dimension + j],2);
                //e += std::pow(d.f(grid[i * grid_dimension + j]) - u1[i * grid_dimension + j], 2);
            }
        }


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
        //MPI_Allreduce(MPI_IN_PLACE, &e, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    // il seguente pezzo di codice serve solo per misurare il tempo impiegato dall'esecuzione
    if (par)
        MPI_Barrier(MPI_COMM_WORLD);
    
    double t_end = MPI_Wtime();

    double local_time = t_end - t_start;

    if (par)
        MPI_Allreduce(MPI_IN_PLACE, &local_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    //if (rank == 0) {
    //    std::cout << "Execution time: " << local_time << " seconds\n";
    //}

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
            int idx = i * n + j; // Adatta questo indice alla struttura della tua matrice
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




/* std::vector<double> solver::solver< s, m, bc, proc, thr>::set_boundary_conditions( ... ) {
    if (constexpr m == DIRICHLET) {
        u = set_dirichlet_boundary_conditions( ... )
    }
    if (constexpr m == NEUMANN) {
        u = set_neumann_boundary_conditions( ... )
    }
    if (constexpr m == ROBIN) {
        u = set_robin_boundary_conditions( ... )
    }

}

std::vector<double> solver::solver< s, m, bc, proc, thr>::set_dirichlet_boundary_conditions( ... ) {
    if (rank == 0) {
        for (int i = 0; i < grid_dimension - 1; ++i)
            u[i] = ...
    }
    for (int i = skip_start_row; i < local_rows - skip_end_row; ++i) {
        u[grid_dimension * (local_previous_rows + i)] = ...
        u[grid_dimension * (local_previous_rows + i) + grid_dimension - 1] = ...
    }
    if (rank == size - 1) {
        for (int i = 0; i < grid_dimension - 1; ++i)
            u[grid_dimension * (grid_dimension - 1) + i] = ...
    }
}

std::vector<double> solver::solver< s, m, bc, proc, thr>::set_neumann_boundary_conditions( ... ) {
    if (rank == 0) {
        u1[0] = 0.25 * (2 * u0[1] + 2 * u0[grid_dimension] - 2 * h * (d.g(grid[ ... ]) + d.g(grid[ ... ]))) - h * h * d.f(grid[ ... ]));

        for (int i = 1; i < grid_dimension - 1; ++i) 
            u1[i] = 0.25 * (2 * u0[i + grid_dimension] - 2 * h * d.g(grid[ ... ]) + u0[i - 1] + u0[i + 1] - h * h * d.f(grid[ ... ]));

        u1[grid_dimension - 1] = 0.25 * (2 * u0[grid_dimension - 2] + 2 * u0[grid_dimension + grid_dimension - 1] + 2 * h * (d.g(grid[ ... ]) - d.g(grid[ ... ]))) - h * h * d.f(grid[ ... ]));
    }
    
    for (int i = skip_start_row; i < local_rows - skip_end_row; ++i) {
        u1[grid_dimension * (local_previous_rows + i)] = 0.25 * (2 * u0[grid_dimension * (local_previous_rows + i) + 1] - 2 * h * d.g(grid[ ... ]) + u0[grid_dimension * (local_previous_rows + i - 1)] + u0[grid_dimension * (local_previous_rows + i + 1)] - h * h * d.f(grid[ ... ]));
        u1[grid_dimension * (local_previous_rows + i) + grid_dimension - 1] = 0.25 * (2 * u0[grid_dimension * (local_previous_rows + i) + grid_dimension - 2] + 2 * h * d.g(grid[ ... ]) + u0[grid_dimension * (local_previous_rows + i - 1) + grid_dimension - 1] + u0[grid_dimension * (local_previous_rows + i + 1) + grid_dimension - 1] - h * h * d.f(grid[ ... ]));
    }

    if (rank == size - 1) {
        u1[grid_dimension * (grid_dimension - 1)] = 0.25 * (2 * u0[grid_dimension * (grid_dimension - 1) + 1] + 2 * u0[grid_dimension * (grid_dimension - 2)] - 2 * h * (d.g(grid[ ... ]) - d.g(grid[ ... ]))) - h * h * d.f(grid[ ... ]));

        for (int i = 1; i < grid_dimension - 1; ++i)
            u1[grid_dimension * (grid_dimension - 1) + i] = 0.25 * (2 * u0[grid_dimension * (grid_dimension - 1) + i - grid_dimension] + u0[grid_dimension * (grid_dimension - 1) + i - 1] + u0[grid_dimension * (grid_dimension - 1) + i + 1] + 2 * h * d.g(grid[ ... ]) - h * h * d.f(grid[ ... ]));

        u1[grid_dimension * grid_dimension - 1] = 0.25 * (2 * u0[grid_dimension * grid_dimension - 2] + 2 * u0[grid_dimension * (grid_dimension - 1) - grid_dimension] + 2 * h * (d.g(grid[ ... ]) + d.g(grid[ ... ]))) - h * h * d.f(grid[ ... ]));
    }
}

*/