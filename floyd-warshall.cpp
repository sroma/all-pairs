// Floyd - Warshall algorithm 
// Parallel with OpenMP
// Graph data stores as connectivity matrix (A)
// Infinity = INT_MAX
#include <iostream>
#include <fstream>
#include <climits>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>

using namespace std;

struct Graph {
    int ** A;
    int VortexCount;
    int EdgeCount;
};

// Floyd-Warshall algorithm
// D - Matrix of Distances;
// P - Predessesor Matrix;
// N - Number of vertexes.
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

// Construct path between i and j vertexes
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

    int **P, **D;
    Graph G;
    bool declared = false;
    string GraphType, token;
    vector<int> path; 

    cout << endl << "Roma graph" << endl;
    ifstream fl("./roma99.gr");
    for(string line; getline(fl, line);) {
        if (line[0] == 'c') continue;
        if (line[0] == 'p') {
            istringstream iss(line);
            iss >> token >> GraphType >> G.VortexCount >> G.EdgeCount;

            if (G.VortexCount > 0) {
                D = new int*[G.VortexCount];
                P = new int*[G.VortexCount];
                G.A = new int*[G.VortexCount];
                for (int i = 0; i < G.VortexCount; i++) {
                    D[i] = new int[G.VortexCount];
                    P[i] = new int[G.VortexCount];
                    G.A[i] = new int[G.VortexCount];
                }
                #pragma omp parallel for
                for (int i = 0; i < G.VortexCount; i++)
                    for (int j = 0; j < G.VortexCount; j++) G.A[i][j] = 0;
                declared = true;
            }
        } 
        if (declared && line[0] == 'a') {
            int i, j, w;
            istringstream iss(line);
            iss >> token >> i >> j >> w; 
            G.A[i-1][j-1] = w;     
        }
    }    
    if (!declared) return 1;
     // initialization of Predessesor and Distance Matrix
     # pragma omp parallel for
     for (int i = 0; i < G.VortexCount; i++)
        for (int j = 0; j < G.VortexCount; j++) {
            D[i][j] = G.A[i][j];
            P[i][j] = i;
            if (i != j && D[i][j] == 0) {
                P[i][j] = -INT_MAX;
                D[i][j] = INT_MAX;
            }
        }
    double tm = omp_get_wtime(); 
    FloydWarshall(D, P, G.VortexCount);
    tm = omp_get_wtime() - tm; 
    cout << num_threads << " threads: " << tm << "s for Rome graph." << endl;

    // test route
    cout << "Path 1-788, d = " << D[0][787] << endl;
    ConstructPath(P,0,787,path);    
    for (int i = 0; i < path.size(); i++) cout << path[i] + 1 << "-";
    cout << endl;

    for (int i = 0; i < G.VortexCount; i++) {
        delete[] D[i];
        delete[] P[i];
        delete[] G.A[i];
    }    
    delete[] D; delete[] P; delete[] G.A;
    return 0;
}
