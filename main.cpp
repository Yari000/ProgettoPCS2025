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
        double pathLength=0.0;
        unsigned numedgesOnPath=0;
        if (startv >= mesh.NumVertices || endv >= mesh.NumVertices) {
        throw std::runtime_error("gli indici dei vertici sono fuori dal range");
        }

        computeShortestPath(mesh, startv, endv, Vertonpath, Edgeonpath, pathLength, numedgesOnPath);

        std::cout << "Geodesic mesh generato:\n";
        std::cout << " - Vertici: " << mesh.NumCell0Ds << "\n";
        std::cout << " - Facce triangolari: " << mesh.NumCell2Ds << "\n";

        std::cout << "\nEsempio primo vertice: "
            << mesh.Cell0DsCoordinates.row(0) << "\n";
        std::cout <<"Lunghezza del percorso più breve: "<<pathLength<<"\n";
 


        {
            std::ofstream vout ("vertices.csv", std::ofstream::out);
            vout << "id,x,y,z,sp\n";
            for (unsigned i = 0; i < mesh.NumVertices; ++i) {
                auto p0 = mesh.Cell0DsCoordinates.row(i);
                vout << i << "," << p0.x() << "," << p0.y() << "," << p0.z() << "," << (Vertonpath.[i] ? "true" : "false") << "\n";
            }
        }

        std::ofstream eout ("edges.csv", std::ofstream::out);
        eout << "id,v0,v1,sp\n";
        for (unsigned int i = 0; i < mesh.NumEdges; ++i) {
            unsigned a = mesh.Cell1DsExtrema(i, 0);
            unsigned b = mesh.Cell1DsExtrema(i, 1);
            eout << i << "," << a << "," << b << "," << (Edgeonpath[i] ? "true" : "false") << "\n";
        }

      std::ofstream fout("Cell2Ds.csv", std::ofstream::out);
fout << "id,vertices,edges\n";
for (unsigned int i = 0; i < mesh.NumFaces; ++i) {
    fout << i << ",";
    for (unsigned j = 0; j < mesh.NumVerticesPerCell2D[i]; ++j) {
        fout << mesh.Cell2DsVertices[i][j];
        if (j < mesh.NumVerticesPerCell2D[i] - 1) {
            fout << " ";
        }
    }
	fout << ",";
    for (unsigned j = 0; j < mesh.NumEdgesPerCell2D[i]; ++j) {
        fout << mesh.Cell2DsEdges[i][j];
        if (j < mesh.NumEdgesPerCell2D[i] - 1) {
            fout << " ";
		}
    }
    fout << "\n";
}

std::ofstream meshout("Cell3Ds.csv", std::ofstream::out);
meshout << "id,NumVertices,vertices,NumEdges,edges,NumFaces,faces\n";
for (unsigned int i = 0; i < mesh.NumCell3Ds; ++i) {
    meshout << i << ",";
    meshout << mesh.NumVertices << ",";
    for (unsigned j = 0; j < mesh.Cell3DsVertices[i].size(); ++j) {
        meshout << mesh.Cell3DsVertices[i][j];
        if (j < mesh.Cell3DsVertices[i].size() - 1) {
            meshout << " ";
        }
    }
    meshout << ",";
    meshout << mesh.NumEdges << ",";
    for (unsigned j = 0; j < mesh.Cell3DsEdges[i].size(); ++j) {
        meshout << mesh.Cell3DsEdges[i][j];
        if (j < mesh.Cell3DsEdges[i].size() - 1) {
            meshout << " ";
        }
    }
	meshout << ",";
    meshout << mesh.NumFaces << ",";
    for (unsigned j = 0; j < mesh.Cell3DsFaces[i].size(); ++j) {
        meshout << mesh.Cell3DsFaces[i][j];
        if (j < mesh.Cell3DsFaces[i].size() - 1) {
            meshout << " ";
        }
    }
    meshout << "\n";

}

   Gedim::UCDUtilities utilities;

{
    vector<Gedim::UCDProperty<double>> Cell0Ds_properties(1);
    Cell0Ds_properties[0].Label ="ShortestPath";
    Cell0Ds_properties[0].UnitLabel ="-";
    Cell0Ds_properties[0] =NumComponents =1;

    vector<double> shortest_path_marker(mesh.NumVertices,0);
    for (unsigned i = 0; i < mesh.NumVertices; ++i)
    shortest_path_marker[i] = Vertonpath[i] ? 1.0 : 0.0;

    Cell0Ds_properties[0].Data =shortest_path_marker.data();
    Eigen::VectorXi materials = Eigen::VectorXi::Zero(mesh.NumVertices);
    utilities.ExportPoints("./Cell0Ds.inp", mesh.Cell0DsCoordinates.transpose(), Cell0Ds_properties,materials);
}
  
{
    vector<Gedim::UCDProperty<double>> Cell0Ds_properties(1);
    Cell0Ds_properties[0].Label = "ShortestPath";
    Cell0Ds_properties[0].UnitLabel = "-";
    Cell0Ds_properties[0].NumComponents = 1;
    vector<double> vertex_marker(mesh.NumVertices, 0.0);
    for (unsigned i = 0; i < mesh.NumVertices; ++i) {
        vertex_marker[i] = Vertonpath[i] ? 1.0 : 0.0;
    }
    Cell0Ds_properties[0].Data = vertex_marker.data();
    vector<Gedim::UCDProperty<double>> Cell1Ds_properties(1);
    Cell1Ds_properties[0].Label = "ShortestPath";
    Cell1Ds_properties[0].UnitLabel = "-";
    Cell1Ds_properties[0].NumComponents = 1;
    vector<double> edge_marker(mesh.NumEdges, 0.0);
    for (unsigned i = 0; i < mesh.NumEdges; ++i) {
        edge_marker[i] = Edgeonpath[i] ? 1.0 : 0.0;
    }

    Cell1Ds_properties[0].Data = edge_marker.data();
    Eigen::VectorXi materials = Eigen::VectorXi::Zero(mesh.NumEdges);

    utilities.ExportSegments("./Cell1Ds.inp", mesh.Cell0DsCoordinates.transpose(), mesh.Cell1DsExtrema.transpose(), Cell0Ds_properties, Cell1Ds_properties, materials);

}
    

    }
    catch (const std::exception& e) {
        std::cerr << "Errore: " << e.what() << "\n";     
    }

    return 0;
}

