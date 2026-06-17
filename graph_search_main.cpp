// Compilation: g++ -std=c++11 -fopenmp -o graph.exe
// Start ./graph.exe 20
#include <list>
#include <vector>
#include <iostream>
#include <memory>
#include <fstream>
#include <climits>
#include <string>
#include <set>
#include <utility>
#include <sstream>
#include <exception>
#include <omp.h>
using namespace std;

struct Edge {
    Edge(const int vNum, const int eWeight): v(vNum), w(eWeight) {}
    int v, w;
};

typedef list <Edge> List;
typedef vector <List> Lists;

struct Graph {
    Graph(const int vCount, const int eCount):
     VertexCount(vCount), EdgeCount(eCount), pLists(new Lists(vCount)) {}
    int VertexCount, EdgeCount;
    unique_ptr<Lists> pLists;
};

Graph loadFromFile(const char * name) {
    int vCount, eCount;
    string graphType, token;
    ifstream f(name);
    for(string line; getline(f,line);) {
        if (line[0] == 'c') continue;
        if (line[0] == 'p') {
            istringstream iss(line);
            iss >> token >> graphType >> vCount >> eCount;
        }
    };
    if (vCount <= 0 || eCount <= 0) throw invalid_argument("No graph data is found.");

    Graph G = Graph(vCount, eCount);
    f.clear();
    f.seekg(0, ios::beg);
    for(string line; getline(f,line);) {
        if (line[0] == 'c' || line[0] == 'p') continue;
        if (line[0] == 'a') {
            int u, v, w;
            istringstream iss(line);
            iss >> token >> u >> v >> w;
            if(u <= G.VertexCount && v <= G.VertexCount)
             G.pLists->at(u-1).push_back(Edge(v-1, w));
        }
    }
    f.close();
    return G;
}

list<int> ConstructPath(int ** Predessesors, int i, int j) {
    list<int> path;
    if (Predessesors[i][j] == -INT_MAX) {
        path.push_back(-INT_MAX);
        return path;
    }
    while (i != j) {
        path.push_front(j);
        j = Predessesors[i][j];
    }
    path.push_front(i);
    return path;
}

void init(Graph * G, int** Distances, int** Predessesors) {
    # pragma omp parallel for
    for (auto i = 0; i < G->VertexCount; i++) {
        for (auto j = 0; j < G->VertexCount; j++) {
            Distances[i][j] = INT_MAX;
            Predessesors[i][j] = INT_MAX;
        }
        for (auto e: G->pLists->at(i)) {
            Distances[i][e.v] = e.w;
            Predessesors[i][e.v] = i;
        }
        Distances[i][i] = 0;
    }
}

void FloydWarshall(Graph * G, int** Distances, int** Predessesors) {
    init(G, Distances, Predessesors);
    for (int k = 0; k < G->VertexCount; k++) {
        #pragma omp parallel for
        for (int i = 0; i < G->VertexCount; i++)
            for (int j = 0; j < G->VertexCount; j++)
                if (Distances[i][k] < INT_MAX && Distances[k][j] < INT_MAX) {
                    int d = Distances[i][k] + Distances[k][j];
                    if (d < Distances[i][j]) {
                         Distances[i][j] = d;
                         Predessesors[i][j] = Predessesors[k][j];
                    }
                }
    }
}

bool BellmanFord(Graph * G, int Source, int * Distances, int * Predessesors) {
    Distances[Source] = 0;
    for (int iter = 0; iter <= G->VertexCount; iter++) // last iter checks negative cycle
        for (int i = 0; i < G->VertexCount; i++)
            for (auto e: G->pLists->at(i)) {
                int d = Distances[i] + e.w;
                if (Distances[i] < INT_MAX && Distances[e.v] > d) {
                    if (iter == G->VertexCount) return true;
                    Distances[e.v] = d;
                    Predessesors[e.v] = i;
                }
            }
    return false;
}

// Dijkstra alg. (array version)
void Dijkstra_arr(Graph * G, int Source, int * Distances, int * Predessesors) {
    bool * u = new bool[G->VertexCount];
    #pragma omp parallel for
    for (int i = 0; i < G->VertexCount; i++) u[i] = false;
    Distances[Source] = 0;

    for (int i = 0; i < G->VertexCount; i++) {
        int mind = INT_MAX;
        int cur = INT_MAX;
        for (int j = 0; j < G->VertexCount; j++)
            if (Distances[j] < mind && !u[j]) {
                mind = Distances[j];
                cur = j;
            }

        if (mind == INT_MAX || cur == INT_MAX) break; // graph has non-reacheable vortex

        u[cur] = true;
        for (auto e: G->pLists->at(cur)) {
                int d = Distances[cur] + e.w;
                if (Distances[e.v] > d) {
                    Distances[e.v] = d;
                    Predessesors[e.v] = cur;
                }
            }
    }
}
// Dijkstra alg. (associate arrays: get perfomance of O(VlogE) instead of O(V^2))
void Dijkstra_set(Graph * G, int Source, int * Distances, int * Predessesors) {
    set<pair<int, int>> unvisited; // set of unvisited vortices
    // start at Source
    Distances[Source] = 0;
    // none is reachable from Source
    for (auto e: G->pLists->at(Source)) Distances[e.v] = INT_MAX;
    unvisited.emplace(0, Source);

     while (!unvisited.empty()) {
         // get vortex with min distance (it is the first in set)
         int cur = unvisited.begin()->second;
         unvisited.erase(unvisited.begin());
         // find all neighbors of cur
         for (auto e: G->pLists->at(cur)) {
            int d = Distances[cur] + e.w;
            if (d < Distances[e.v]) {
                // remove visited vortex
                unvisited.erase({Distances[e.v], e.v});
                // update distance
                Distances[e.v] = d;
                Predessesors[e.v] = cur;
                unvisited.emplace(d,e.v);
            }
        }
    } // end while
}

// Johnson alg. with Dijkstra alg. using arrays
bool Johnson(Graph * G, int** Distances, int** Predessesors) {
    init(G, Distances, Predessesors);

    if (BellmanFord(G, 0, Distances[0], Predessesors[0])) return true;
    #pragma omp parallel for
    for (int i = 0; i < G->VertexCount; i++)
        Dijkstra_arr(G, i, Distances[i], Predessesors[i]);
    return false;
}

// Johnson alg. with Dijkstra alg. using associate arrays
bool Johnson2(Graph * G, int** Distances, int** Predessesors) {
    init(G, Distances, Predessesors);

    if (BellmanFord(G, 0, Distances[0], Predessesors[0])) return true;
    #pragma omp parallel for
    for (int i = 0; i < G->VertexCount; i++)
        Dijkstra_set(G, i, Distances[i], Predessesors[i]);
    return false;
}

// test case
int main(int argc, char* argv[]) {
    int num_threads = (argc > 1) ? atoi(argv[1]) : 1;
    omp_set_num_threads(num_threads);

    Graph G = loadFromFile("1.txt");
    int** Distances = new int*[G.VertexCount];
    int** Predessesors = new int*[G.VertexCount];
    for (int i = 0; i < G.VertexCount; i++) {
        Distances[i] = new int[G.VertexCount];
        Predessesors[i] = new int[G.VertexCount];
    }
    cout << "Threads: " << num_threads << "; V = " <<
        G.VertexCount << "; E = " << G.EdgeCount << endl;

    // Floyd-Warshall test
    double tm = omp_get_wtime();
    FloydWarshall(&G, Distances, Predessesors);
    cout << "Floyd-Warshall: "<< omp_get_wtime() - tm << "s." << endl;
    cout << "Distance 955-315: " << Distances[954][314] << endl;
    list<int> path = ConstructPath(Predessesors,954,314);
    for (auto v: path) cout << v + 1 << "-";
    cout << endl;

    // Johnson + Dijkstra (associate arr.)
    tm = omp_get_wtime();
    if (Johnson2(&G, Distances, Predessesors)) {
        cout << "Negative cycle is in graph" << endl;
    }
    else {
        cout << "Johnson (set in Dijkstra): "<< omp_get_wtime() - tm << "s." << endl;
        cout << "Distance 955-315: " << Distances[954][314] << endl;
        path = ConstructPath(Predessesors,954,314);
        for (auto v: path) cout << v + 1 << "-";
        cout << endl;
    }

    // Johnson + Dijkstra (arr.)
    tm = omp_get_wtime();
    if (Johnson(&G, Distances, Predessesors)) {
        cout << "Negative cycle in graph" << endl;
    }
    else {
        cout << "Johnson (array in Dijkstra): "<< omp_get_wtime() - tm << "s." << endl;
        cout << "Distance 955-315: " << Distances[954][314] << endl;
        path = ConstructPath(Predessesors,954,314);
        for (auto v: path) cout << v + 1 << "-";
        cout << endl;
    }

    // cleaning
    for (int i = 0; i < G.VertexCount; i++) {
        delete[] Distances[i];
        delete[] Predessesors[i];
    }
    delete[] Distances;
    delete[] Predessesors;
  
    return 0;
}
