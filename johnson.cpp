#include <iostream>
#include <fstream>
#include <climits>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <omp.h>
using namespace std;

struct Edge {
    int start;
    int end;
    int weight;
};
struct Graph {
    Edge * E;
    int VortexCount;
    int EdgeCount;
};

void ConstructPath(int j, int * p, vector<int>& path) {
    path.clear();
    if (p[j] == INT_MAX) path.push_back(-INT_MAX);
    else 
        for (int i = j; i != -1; i = p[i])
            path.push_back(i);
}

void PrintPath(Graph G, int i, int j, int d, int * p) {
    vector<int> path;
    cout << "Path " << i+1 << "-" << j+1 << ", d = " << d << endl;
    ConstructPath(j, p, path);    
    reverse(path.begin(), path.end());
    for (int i = 0; i < path.size(); i++) cout << path[i] + 1 << "-";
    cout << endl;
    int s = 0;
    for (int i = 0; i < path.size() - 1; i++)
	for (int j = 0; j < G.EdgeCount; j++) 
	    if (G.E[j].start == path[i] && G.E[j].end == path[i+1]) {
    		cout << G.E[j].weight << "+";
    		s +=  G.E[j].weight;
    		break;
    	    }
    cout << " = " << s << endl;
}

// Bellman-Ford algorithm, s - source of search
// d - length from k to other vertexes
bool BellmanFord(Graph G, int s, int * d, int * p) {
    #pragma omp parallel for
    for (int i = 0; i < G.VortexCount; i++) {
        d[i] = INT_MAX;
        p[i] = -1;
    }
    d[s] = 0;
    for (int i = 0; i <= G.VortexCount; i++) // last iter check negative cycle
	    for (int j = 0; j < G.EdgeCount; j++) {
	        if (d[G.E[j].start] < INT_MAX && d[G.E[j].end] > d[G.E[j].start] + G.E[j].weight) {
		        if (i == G.VortexCount) return true;
                d[G.E[j].end] = d[G.E[j].start] + G.E[j].weight;
                p[G.E[j].end] = G.E[j].start;
            }    
        }
    return false;
}

// Dijkstra alg (marks and distances store in arrays)
void Dijkstra(Graph G, int s, int *d, int * p) {
    bool * u = new bool[G.VortexCount];
    #pragma omp parallel for
    for (int i = 0; i < G.VortexCount; i++) {
        d[i] = INT_MAX;
        p[i] = -1;
	    u[i] = false;
    }    
    d[s] = 0;
    
    for (int i = 0; i < G.VortexCount; i++) {
	    int mind = INT_MAX;
	    int cur = INT_MAX;
	    for (int j = 0; j < G.VortexCount; j++) 
	        if (d[j] < mind && !u[j]) {
		        mind = d[j];
		        cur = j;
	        }
	
	    if (mind == INT_MAX || cur == INT_MAX) break; // graph is not connected
	    
	    u[cur] = true;
	    for (int j = 0; j < G.EdgeCount; j++) {
	        if (G.E[j].start == cur) {
		        int dist = d[cur] + G.E[j].weight;
		        if (d[G.E[j].end] > dist) {
		            d[G.E[j].end] = dist;
		            p[G.E[j].end] = cur;
		        }
	        }
	    }
    }
}

int main(int argc, char* argv[]) {
    int num_threads = (argc > 1) ? atoi(argv[1]) : 1; 
    omp_set_num_threads(num_threads);

    Graph G;
    int num = 0;
    string GraphType, token;
    bool declared = false;
    cout << endl << "Roma graph" << endl;
    ifstream fl("./roma99.gr");
    for(string line; getline(fl, line);) {
        if (line[0] == 'c') continue;
        if (line[0] == 'p') {
            istringstream iss(line);
            iss >> token >> GraphType >> G.VortexCount >> G.EdgeCount;
	    
	        if (G.VortexCount > 0 && G.EdgeCount > 0) {
		        G.E = new Edge[G.EdgeCount];
	    	    declared = true;
	        }
	    
        } 
        if (declared && line[0] == 'a') {
            int i, j, w;
            istringstream iss(line);
            iss >> token >> i >> j >> w; 
            G.E[num].start = i-1;
	        G.E[num].end = j-1;
	        G.E[num].weight = w;
	        num++;     
        }
    } 
          
    double tm = omp_get_wtime(); 

    int ** D = new int*[G.VortexCount + 1];
    int ** P = new int*[G.VortexCount + 1];
    for (int i = 0; i < G.VortexCount + 1; i++) {
	    D[i] = new int [G.VortexCount + 1];
	    P[i] = new int [G.VortexCount + 1];
    }
    Graph G2 = G;
    G2.VortexCount = G.VortexCount + 1;
    G2.EdgeCount = G.EdgeCount + G.VortexCount;
    G2.E = new Edge[G2.EdgeCount];
    
    #pragma omp parallel for
    for (int i = 0; i < G.VortexCount; i++) G2.E[i] = G.E[i];
    #pragma omp parallel for
    for (int i = 0; i < G2.VortexCount-1; i++) {
	    G2.E[G.EdgeCount + i].start = G.VortexCount;
	    G2.E[G.EdgeCount + i].end = i;
	    G2.E[G.EdgeCount + i].weight = 0;
    }	
    
    BellmanFord(G2, G.VortexCount, D[G.VortexCount], P[G.VortexCount]); 
    #pragma omp parallel for
    for (int i = 0; i < G.EdgeCount; i++) G.E[i].weight -= (D[G.VortexCount][G.E[i].start] - D[G.VortexCount][G.E[i].end]);
    #pragma omp parallel for
    for (int i = 0; i < G.VortexCount; i++) Dijkstra(G,i,D[i],P[i]);
    #pragma omp parallel for
    for (int i = 0; i < G.EdgeCount; i++) G.E[i].weight += (D[G.VortexCount][G.E[i].start] - D[G.VortexCount][G.E[i].end]);
    
    cout << num_threads << " threads: " << omp_get_wtime() - tm << "s for Rome graph." << endl;
    
    // Test case
    //for (int i = 0; i < G.VortexCount; i++)
    //   cout << 0 << "-" << i << " " << D[0][i] - (D[G.VortexCount][0 - D[i]) << endl;
    PrintPath(G, 0, 787, D[0][787] - (D[G.VortexCount][0 - D[787]) , P[0]); 
    PrintPath(G, 0, 147, D[0][147] - (D[G.VortexCount][0 - D[147]) , P[0]); 

    delete[] G.E;
    delete[] G2.E;
    for (int i = 0; i < G.VortexCount + 1; i++) {
	    delete[] D[i];
	    delete[] P[i];
    }
    delete[] D;
    delete[] P;
    return 0;
}
