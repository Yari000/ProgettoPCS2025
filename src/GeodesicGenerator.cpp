#include "GeodesicGenerator.hpp"
#include "GeometryUtils.hpp"
#include "PolygonalMesh.hpp"
#include <map>
#include <unordered_map>
#include <fstream>
#include <queue>
#include <set>
#include <unordered_set>
#include <stdexcept>


namespace PolygonalLibrary {

    PolygonalMesh generateGeodesicMesh(unsigned int p, unsigned int q, unsigned int b) {

        // genera la mesh basata sui solidi platonici di base
        PolygonalMesh mesh;
        if (p != 3 and q != 3) {
            throw std::runtime_error("solo solidi con facce triangolari e loro duali sono supportati (p,q=3)");
        }

        std::vector<Eigen::Vector3d> verts;
        std::vector<std::array<unsigned, 3>> tris;
		std::vector<Eigen::Vector3d> overts; // vertici della mesh dopo la suddivisione
		std::vector<std::vector<unsigned>> otris; // facce della mesh dopo la suddivisione


        if (q==5 or p==5) {
            generateIcosahedron(verts, tris);
        }
        else if (q==3 and p==3) {
            generateTetrahedron(verts, tris);
        }
        else if (q==4 or p==4) {
            generateOctahedron(verts, tris);
        }
        else {
            throw std::runtime_error("p,q deve essere 3, 4 o 5 per tetraedro, ottaedro o icosaedro e relativi duali");
        }



        // suddivisione della mesh
        subdivideGeometry_T(verts, tris, b,overts,otris);

        
        // calcolo del duale della mesh
        if (q == 3 and p != 3) {
            std::vector<Eigen::Vector3d> dualVerts;
            std::vector<std::vector<unsigned>> dualFaces;
            computeDualMesh(overts, otris, dualVerts, dualFaces);

            overts = dualVerts;
	    otris = dualFaces;
            }
        }
        
        if (p == 3) {
            std::vector<std::vector<unsigned>> polyFaces;
            for (const auto& tri : otris) {
                polyFaces.push_back({ tri[0], tri[1], tri[2] });
            }
            otris = polyFaces;
        }
        
        

        // assegnazione dei dati alla structure della mesh

        //Cell0Ds
        mesh.NumCell0Ds = overts.size();
        mesh.Cell0DsCoordinates.resize(mesh.NumCell0Ds, 3);
        for (unsigned i = 0; i < overts.size(); ++i) {
            mesh.Cell0DsCoordinates.row(i) = overts[i].transpose();
            mesh.Cell0DsId.push_back(i);
        }
        
        //Cell2Ds
        mesh.NumCell2Ds = otris.size();
        mesh.Cell2DsId.resize(otris.size());
        mesh.NumVerticesPerCell2D.resize(otris.size());
	mesh.NumEdgesPerCell2D.resize(otris.size());
        mesh.Cell2DsVertices.resize(otris.size());
        for (unsigned i = 0; i < otris.size(); ++i) {
            const auto& face = otris[i];
            mesh.Cell2DsId[i] = i;
            mesh.NumVerticesPerCell2D[i] = face.size();
            mesh.NumEdgesPerCell2D[i] = face.size();
            mesh.Cell2DsVertices[i] = face;
        }
            
            
        //Cell3Ds
        mesh.NumVertices = mesh.NumCell0Ds;
        mesh.NumFaces = mesh.NumCell2Ds;

        mesh.NumCell3Ds= 1;
        mesh.Cell3DsId= {0};
        mesh.Cell3DsVertices.resize(1);
	mesh.Cell3DsEdges.resize(1);
	mesh.Cell3DsFaces.resize(1);

        for (unsigned i = 0; i < mesh.NumCell0Ds; ++i) {
            mesh.Cell3DsVertices[0].push_back(i);
		}
        
        // Cell1Ds
        std::map<std::pair<unsigned, unsigned>, unsigned> edgeMap;
        mesh.Cell1DsId.clear();
	mesh.Cell1DsExtrema.resize(0, 2);

	unsigned edgeCount= 0;
        for (unsigned f = 0; f < otris.size(); ++f) {
            const auto& face = otris[f];
            std::vector<unsigned> edgeIds;
            for (unsigned i = 0; i < face.size(); ++i) {
                unsigned a = face[i];
                unsigned b = face[(i + 1) % face.size()];
                if (a > b) std::swap(a, b); // ordina gli estremi dell'edge per evitare duplicati

                auto edge = std::make_pair(a, b);
                if (edgeMap.count(edge) == 0) {
                    edgeMap[edge] = edgeCount++;
                    mesh.Cell1DsId.push_back(edgeMap[edge]);
                }
                edgeIds.push_back(edgeMap[edge]);
            }
            mesh.Cell2DsEdges.push_back(edgeIds);
        }

        // costruzione degli estremi dei lati
		mesh.NumCell1Ds = edgeMap.size();
		mesh.Cell1DsExtrema.resize(mesh.NumCell1Ds, 2);
        for (const auto& [edge, id] : edgeMap) {
            mesh.Cell1DsExtrema(id, 0) = edge.first;
            mesh.Cell1DsExtrema(id, 1) = edge.second;
		}

		mesh.NumEdges = mesh.NumCell1Ds;

		
                for (unsigned i = 0; i < mesh.NumCell1Ds; ++i) {
			mesh.Cell3DsEdges[0].push_back(i);
		}
		for (unsigned i = 0; i < mesh.NumCell2Ds; ++i) {
			mesh.Cell3DsFaces[0].push_back(i);  
		}


        return mesh;

    }



}

