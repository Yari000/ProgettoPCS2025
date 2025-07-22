#include "ShortestPath.hpp"
#include "PolygonalMesh.hpp"
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <utility>


namespace std{
    template <>
    struct hash<std::pair<int, int>> {
        std::size_t operator()(const std::pair<int, int>& p) const {
			std::size_t h1 = std::hash<int>()(p.first);
			std::size_t h2 = std::hash<int>()(p.second);
			return h1 ^ (h2 << 1); // combina gli hash dei due elementi
        }
    };
}


using namespace PolygonalLibrary;



double edgeLength(const Eigen::MatrixXd& coords, int v1, int v2) {
        Eigen::Vector3d p1 = coords.row(v1);
        Eigen::Vector3d p2 = coords.row(v2);
        return (p1 - p2).norm();
}



void computeShortestPath(
    PolygonalMesh& mesh,
	unsigned int startVertex,
	unsigned int endVertex,
	std::vector<bool>& vertexOnPath,
	std::vector<bool>& edgeOnPath)

{ 

const auto& coords = mesh.Cell0DsCoordinates; 

// crea una matrice di adiacenza per il grafo
std::unordered_map<int, std::vector<std::pair<int, double>>> graph;
for (int i = 0; i < mesh.NumCell1Ds; ++i) {
    unsigned v1 = mesh.cell1DsExtrema(i, 0);
    unsigned v2 = mesh.cell1dsExtrema(i, 1);
	double length = edgeLength(coords, v1, v2);
    graph[v1].emplace_back(v2, length);
    graph[v2].emplace_back(v1, length);   // aggiunge la distanza  tra i vertici alla matrice
}


// implementazione dell'algoritmo di Dijkstra per trovare il cammino più breve
std::vector<double> dist(mesh.NumVertices, std::numeric_limits<double>::infinity());
std::vector<int> prev(mesh.NumVertices, -1);
dist[startVertex] = 0.0;

std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> queue;
queue.push({ 0.0, startVertex }); // inizia da un vertice di partenza

while (!queue.empty()) {
    auto [d, u] = queue.top(); queue.pop();

	if (u == endVertex) break; // se raggiunge il vertice di arrivo, esce
    if (d > dist[u]) continue; // se la distanza è già maggiore, salta

    for (auto [v, len] : graph[u]) {
        if (dist[u] + len < dist[v]) {
            dist[v] = dist[u] + len;
            prev[v] = u;
            queue.push({ dist[v], v });
        }
    }
}


// trova il cammino più breve tra il vertice di partenza e quello di arrivo
vertexOnPath.resize(mesh.NumVertices, false);
edgeOnPath.resize(mesh.NumEdges, false);

for (int v=endVertex; v != -1; v = prev[v]) {
    vertexOnPath[v] = true; // segna i vertici sul cammino
}

for (int v = endVertex; prev[v] != -1; v = prev[v]) {
	int u = prev[v]; // il vertice precedente
    
    for (int e=0; e<mesh.NumEdges; ++e) {
            int a = mesh.Cell1DsExtrema(e, 0);
            int b = mesh.Cell1DsExtrema(e, 1);
            if ((a == u && b == v) || (a == v && b == u)) {
                edgeOnPath[e] = true; // segna gli edge sul cammino
                break;
			}
        }
	}
}