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

// Floyd-Warshall
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

// get path
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

void ConstructPath6(int P[6][6], int i, int j) {
    if (i == j) cout << i+1 << " ";
    else if (P[i][j] == -INT_MAX)
        cout << i+1 << " - " << j+1;
    else {
        ConstructPath6(P, i, P[i][j]); 
        cout << j+1 << " ";
    }
}

int main(int argc, char* argv[]) {
    /* int A6[6][6] = {{0,10,20,30,0,0},
                      {0,0,0,0,0,7},
                      {0,0,0,0,0,5},
                      {0,0,0,0,10,0},
                      {2,0,0,0,0,4},
                      {0,5,7,0,6,0}};
    int D6[6][6], P6[6][6];
     cout << "Test case" << endl;
    // Floyd-Warshall (test case)
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            D6[i][j] = A6[i][j];
            P6[i][j] = i;
            if (i != j && D6[i][j] == 0) {
                D6[i][j] = INT_MAX;
                P6[i][j] = -INT_MAX;
            }
        }
    for (int k = 0; k < 6; k++)
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) 
                if (D6[i][k] < INT_MAX && D6[k][j] < INT_MAX) {
                    int d = D6[i][k] + D6[k][j];
                    if (d < D6[i][j]) {
                         D6[i][j] = d;
                         P6[i][j] = P6[k][j];
                    }
                }

    cout << "Distance matrix" << endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) 
            if (D6[i][j] < INT_MAX) cout << D6[i][j] << " ";
            else cout << "inf ";
        cout << endl;
    }
    cout << "Predessesor matrix" << endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) 
            cout << P6[i][j] << " ";
        cout << endl;
    }
    cout << "Path 1-5, d = " << D6[0][4] << ": ";
    ConstructPath6(P6,0,4); 
    */
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
        } // p
        if (declared && line[0] == 'a') {
            int i, j, w;
            istringstream iss(line);
            iss >> token >> i >> j >> w; 
            G.A[i-1][j-1] = w;     
        }
    }    
    if (!declared) return 1;

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

    //for (int i = 0; i < G.VortexCount; i++)
    //   cout << 0 << "-" << i << " " << D[0][i] << endl;

    cout << "Path 1-788, d = " << D[0][787] << endl;
    ConstructPath(P,0,787,path);    
    for (int i = 0; i < path.size(); i++) cout << path[i] + 1 << "-";
    cout << endl;
    int s = 0;
    for (int i = 0; i < path.size() - 1; i++){
        cout << G.A[path[i]][path[i+1]] << "+";
        s += G.A[path[i]][path[i+1]];
    }
    cout << " = " << s << endl;
    /*
    cout << "Path 1-148, d = " << D[0][147] << endl;
    ConstructPath(P,0,147,path);    
    for (int i = 0; i < path.size(); i++) cout << path[i] + 1 << "-";
    cout << endl;
    s = 0;
    for (int i = 0; i < path.size() - 1; i++){
        cout << G.A[path[i]][path[i+1]] << "+";
        s += G.A[path[i]][path[i+1]];
    }
    cout << " = " << s << endl;
    */
    for (int i = 0; i < G.VortexCount; i++) {
        delete[] D[i];
        delete[] P[i];
        delete[] G.A[i];
    }    
    delete[] D; delete[] P; delete[] G.A;
    return 0;
}
