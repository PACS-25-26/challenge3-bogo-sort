#include "solver.hpp"


#include <functional>
#include <cmath>
#include <vector>
#include <iostream>


Solution solve_pde(ProblemData d){
    // MPI_Init(int argc, char* argv[]);
    // int rank, size;
    // MPI_Comm_rank(&rank);
    // MPI_Comm_size(&size);

    // int n;
    // if (rank == 0) {
    //      std::cout << "Type the number of elements desired: "
    //      std::cin >> n;
    // }
    // sostituire questa riga che permette solo dimensioni pari della griglia
    int grid_dimension = std::pow(2,d.nref)+1;
    // int grid_dimension = n;
    double h = (d.x1-d.x0) / (grid_dimension-1);

    // questi saranno i vettori global
    Vector u1(grid_dimension * grid_dimension,0);
    Vector u0(grid_dimension * grid_dimension,0);
    
    std::cout << grid_dimension << std::endl;
    int iterations = 0;

    // questo sarà l'errore globale
    double error = d.tol + 1;

    // questa sarà la griglia globale
    std::vector<std::vector<double>> grid(grid_dimension * grid_dimension);
    // std::vector<std::vector<double>> global_grid( ... );

    // per implementare la divisione delle righe tra i processi anche se le righe (e quindi le colonne sono dispari), si usa la strategia:
    // int proc = grid_dimension / size;
    // int remainder = grid_dimension % size;
    // std::vector<int> rows_to_rank(size, proc);
    // while (remainder > 0) {
    //      rows_to_rank[remainder] += 1;
    //      remainder--;
    // }
    // MPI_Bcast(&rows_to_rank, 4, MPI_INT, 0, MPI_COMM_WORLD);
    // int local_rows = rows_to_rank[rank];
    //
    // in questo modo creo un vettore contenente quante righe assegnare ad ogni processo e poi broadcasto il vettore. seppur questo vettore non 
    // pesi molto, può essere inutile che tutti i processi conoscano tutte le righe per ogni processo: si può fare solo Send e Receive tipo:
    // int local_rows;
    // for (i = 0; i < size; ++i) {
    //      if (rank == 0) {
    //          MPI_Send(&rows_to_rank[rank - 1], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    //      }
    //      else {
    //          MPI_Recv(&local_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //      }
    // }
    // 
    // int local_cols = grid_dimension;
    // MPI_Bcast(&local_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //
    // 

    for(int i = 0;i < grid_dimension;i++){ // sostituire grid_dimension con local_rows
            for(int j = 0;j < grid_dimension;j++){ // sostituire grid_dimension con local_cols (== grid_dimension)
                grid[i * grid_dimension + j]={d.x0 + i * h, d.y0 + j * h};
            }
        }

    while(error > d.tol and iterations<d.max_iter){
        double sum = 0;
        for(int i = 1;i < grid_dimension - 1;i++){ // for (int i = local_rows * rank + 1; i < local_rows * rank + local_rows - 1; ++i)
            for(int j = 1;j < grid_dimension - 1;j++){ // for (j = 1; j < local_cols; ++j)
                u1[i * grid_dimension + j] = (0.25) * (u0[(i + 1) * grid_dimension + j] +u0[(i - 1) * grid_dimension + j]
                                                             +u0[i * grid_dimension + j + 1] +u0[i * grid_dimension + j - 1] 
                                                             +d.f(grid[i * grid_dimension + j])* h * h);
                sum += std::pow(u1[i * grid_dimension + j] - u0[i * grid_dimension + j],2);
            }
        }
        error = std::sqrt(h*sum);
    //  if (rank != 0) {
    //      MPI_Send(&u1[local_rows * rank * grid_dimension], grid_dimension, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    //      MPI_Recv(&u1[local_rows * rank * grid_dimension - grid_dimension], grid_dimension, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //  }
    // if (rank != size - 1) {
    //      MPI_Send(&u1[local_rows * rank * grid_dimension + (local_rows - 1) * grid_dimension], grid_dimension, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD)
    //      MPI_Recv(&u1[local_rows * rank * grid_dimension + local_rows * grid_dimension], grid_dimension, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
    //  }
        u0=u1;
        iterations++;
    }
    // MPI_Allreduce(&u1, MPI_IN_PLACE, grid_dimension * grid_dimension, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // MPI_Allreduce(&iterations, MPI_IN_PLACE, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    Solution s;
    
    s.u=u1;
    s.n_iter=iterations;
    s.grid=grid;
    s.grid_dimension=grid_dimension;
    s.h=h;
    return s;
    // MPI_Finalize();
}


#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void create_vtk(const std::string& filename, const std::vector<double>& u, int n, double h) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Errore nell'apertura del file VTK!" << std::endl;
        return;
    }

    // 1. Intestazione obbligatoria
    out << "# vtk DataFile Version 3.0\n";
    out << "Soluzione Equazione di Poisson 2D\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";

    // 2. Geometria della griglia (nX, nY, nZ)
    // Usiamo n nodi per X e Y, e 1 solo nodo per Z visto che siamo in 2D
    out << "DIMENSIONS " << n << " " << n << " 1\n";
    out << "ORIGIN 0.0 0.0 0.0\n";
    out << "SPACING " << h << " " << h << " 0.0\n";

    // 3. Dati dei nodi
    out << "POINT_DATA " << n * n << "\n";
    out << "SCALARS u double 1\n";
    out << "LOOKUP_TABLE default\n";

    // 4. Scrittura dei valori
    // ATTENZIONE: VTK vuole prima la X che varia (colonne j) e poi la Y (righe i)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i * n + j; // Adatta questo indice alla struttura della tua matrice
            out << u[idx] << "\n";
        }
    }

    out.close();
    std::cout << "File " << filename << " saved!" << std::endl;
}