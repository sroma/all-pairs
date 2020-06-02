#include <iostream>
#include <fstream>
#include <climits>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>

using namespace std;

struct Graph {
    int ** A; // Connectivetely Matrix
    int VortexCount;
    int EdgeCount;
};

// Floyd-Warshall algorithm
int FloydWarshall(int ** D, int ** P, int N) {
    for (int k = 0; k < N; k++)
        #pragma omp parallel for
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) 
                if (D[i][k] < INT_MAX && D[k][j] < INT_MAX) {
                    int d = D[i][k] + D[k][j];
                    if (d < D[i][j]) {
                         D[i][j] = d;
                         P[i][j] = P[k][j];
                    }
                }
    return 0;
}


void ConstructPath(int ** P, int i, int j, vector<int>& path) {
    if (i == j) {
        path.clear();
        path.push_back(i);
    }
    else if (P[i][j] == -INT_MAX) {
        path.clear();
        path.push_back(-INT_MAX);
        } else {
            ConstructPath(P, i, P[i][j], path); 
            path.push_back(j);
        }
}

int main(int argc, char* argv[]) {
    int num_threads = (argc > 1) ? atoi(argv[1]) : 1; 
    omp_set_num_threads(num_threads);

    int **P, **D; // Predessesors' and Distances' Matrixs
    Graph G;
    bool declared = false;
    string GraphType, token;
    vector<int> path; 


    for (int i = 0; i < G.VortexCount; i++) {
        delete[] D[i];
        delete[] P[i];
        delete[] G.A[i];
    }    
    delete[] D; delete[] P; delete[] G.A;
    return 0;
}
