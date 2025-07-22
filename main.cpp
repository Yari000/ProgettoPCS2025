#include <iostream>
#include "PolygonalMesh.hpp"
#include "GeometryUtils.hpp"
#include "UCDUtilities.hpp"
#include "GeodesicGenerator.hpp"
#include "ShortestPath.hpp"
#include <limits>
#include <cmath>
#include <fstream>
using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;



int main() {
    unsigned int p = 3, q = 4, b = 7, startv=2, endv=5;
    try {
        PolygonalMesh mesh = generateGeodesicMesh(p, q, b);
        std::vector<bool> Vertonpath;
        std::vector<bool> Edgeonpath;
        computeShortestPath(mesh, startv, endv, Vertonpath, Edgeonpath);

        std::cout << "Geodesic mesh generato:\n";
        std::cout << " - Vertici: " << mesh.NumCell0Ds << "\n";
        std::cout << " - Facce triangolari: " << mesh.NumCell2Ds << "\n";

        std::cout << "\nEsempio primo vertice: "
            << mesh.Cell0DsCoordinates.row(0) << "\n";
 


        {
            std::ofstream vout ("vertices.csv", std::ofstream::out);
            vout << "id,x,y,z,sp\n";
            for (unsigned i = 0; i < mesh.NumVertices; ++i) {
                auto p0 = mesh.Cell0DsCoordinates.row(i);
                vout << i << "," << p0.x() << "," << p0.y() << "," << p0.z() << "," << (Vertonpath.count(i) ? "true" : "false") << "\n";
            }
        }

        std::ofstream eout ("edges.csv", std::ofstream::out);
        eout << "id,v0,v1,sp\n";
        for (unsigned int i = 0; i < mesh.NumEdges; ++i) {
            unsigned a = mesh.Cell1DsExtrema(i, 0);
            unsigned b = mesh.Cell1DsExtrema(i, 1);
            auto pr = make_pair(a, b);
            eout << i << "," << a << "," << b << "," << (Edgeonpath.count(pr) ? "true" : "false") << "\n";
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Errore: " << e.what() << "\n";     
    }

    return 0;
}

