#include "solver.hpp"


#include <functional>
#include <cmath>
#include <vector>
#include <iostream>


Solution solve_pde(ProblemData d){
    int grid_dimension = std::pow(2,d.nref)+1;
    double h = (d.x1-d.x0) / (grid_dimension-1);

    
    Vector u1(grid_dimension * grid_dimension,0);
    Vector u0(grid_dimension * grid_dimension,0);
    std::cout << grid_dimension << std::endl;
    int iterations = 0;

    double error = d.tol + 1;

    std::vector<std::vector<double>> grid(grid_dimension * grid_dimension);

    for(int i = 0;i < grid_dimension;i++){
            for(int j = 0;j < grid_dimension;j++){
                grid[i * grid_dimension + j]={d.x0 + i * h, d.y0 + j * h};
            }
        }




    while(error > d.tol and iterations<d.max_iter){
        double sum = 0;
        for(int i = 1;i < grid_dimension - 1;i++){
            for(int j = 1;j < grid_dimension - 1;j++){
                u1[i * grid_dimension + j] = (0.25) * (u0[(i + 1) * grid_dimension + j] +u0[(i - 1) * grid_dimension + j]
                                                             +u0[i * grid_dimension + j + 1] +u0[i * grid_dimension + j - 1] 
                                                             +d.f(grid[i * grid_dimension + j])* h * h);
                sum += std::pow(u1[i * grid_dimension + j] - u0[i * grid_dimension + j],2);
            }
        }
        
        error = std::sqrt(h*sum);

        u0=u1;
        iterations++;
    }
  
    Solution s;
    
    s.u=u1;
    s.n_iter=iterations;
    s.grid=grid;
    s.grid_dimension=grid_dimension;
    s.h=h;
    return s;
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