#include "GeometryUtils.hpp"
#include <map>
#include <cmath>
#include <algorithm>
#include "Eigen/Dense"
#include <limits>


// funzione che genera il tetraedro
void generateTetrahedron(std::vector<Eigen::Vector3d>& verts, std::vector<std::array<unsigned, 3>>& tris) {
    verts.clear();
    tris.clear();

    
    verts.push_back(Eigen::Vector3d(1,1,1));
    verts.push_back(Eigen::Vector3d(-1,-1,1));
    verts.push_back(Eigen::Vector3d(-1,1,-1));
    verts.push_back(Eigen::Vector3d(1,-1,-1));

    tris.push_back({0,1,2});
    tris.push_back({0,3,1});
    tris.push_back({0,2,3});
    tris.push_back({1,3,2});
}

// funzione che genera l'ottedro
void generateOctahedron(std::vector<Eigen::Vector3d>& verts, std::vector<std::array<unsigned, 3>>& tris) {
    verts.clear();
    tris.clear();

    verts.push_back(Eigen::Vector3d(1,0,0));   
    verts.push_back(Eigen::Vector3d(-1,0,0));  
    verts.push_back(Eigen::Vector3d(0,1,0));   
    verts.push_back(Eigen::Vector3d(0,-1,0));  
    verts.push_back(Eigen::Vector3d(0,0,1));   
    verts.push_back(Eigen::Vector3d(0,0,-1));  

    tris.push_back({0,2,4});
    tris.push_back({2,1,4});
    tris.push_back({1,3,4});
    tris.push_back({3,0,4});

    tris.push_back({2,0,5});
    tris.push_back({1,2,5});
    tris.push_back({3,1,5});
    tris.push_back({0,3,5});
}

// funzione che genera l'icosaedro
void generateIcosahedron(std::vector<Eigen::Vector3d>& verts, std::vector<std::array<unsigned, 3>>& tris) {
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    verts = {
        {-1,  phi, 0}, {1,  phi, 0}, {-1, -phi, 0}, {1, -phi, 0},
        {0, -1,  phi}, {0, 1,  phi}, {0, -1, -phi}, {0, 1, -phi},
        { phi, 0, -1}, { phi, 0, 1}, {-phi, 0, -1}, {-phi, 0, 1}
    };
    for (auto& v : verts)
        v.normalize();

    tris = {
        {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
        {1,5,9}, {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
        {3,9,4}, {3,4,2}, {3,2,6}, {3,6,8}, {3,8,9},
        {4,9,5}, {2,4,11}, {6,2,10}, {8,6,7}, {9,8,1}
    };
}

// Suddivide un singolo triangolo in b^2 sotto-triangoli regolari ( suddivisione di tipo 1 )
void subdivideOnce_T(
    const Eigen::Vector3d& A,
    const Eigen::Vector3d& B,
    const Eigen::Vector3d& C,
    int b,
    std::vector<Eigen::Vector3d>& outVerts,
    std::vector<std::vector<unsigned>>& outTris,
    std::map<QuantizedVector, unsigned>& vertexMap,
    double epsilon)
{
    epsilon = 1e-6;
    // Genera i vertici
    std::vector<std::vector<unsigned>> vertGrid(b + 1);

    for (int i = 0; i <= b; ++i) {
        vertGrid[i].resize(i + 1);
        for (int j = 0; j <= i; ++j) {
            double u =1.0 - (double)i / b;
            double v =(double)(i-j) / b;
            double w =(double)j/b;

            Eigen::Vector3d P = u * A + v * B + w * C;
            P.normalize();  // Proietta su sfera unitaria

            // Evita duplicati usando un map
            QuantizedVector qv(P,epsilon);
			auto it = vertexMap.find(qv);
            if (it != vertexMap.end()) {
                vertGrid[i][j] = it->second;
            }
            else {
                unsigned index = outVerts.size();
                outVerts.push_back(P);
                vertexMap[qv] = index;
                vertGrid[i][j] = index;
            }
        }
    }

    // Genera triangoli
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < i; ++j) {
            // Triangolo "sinistro"
            outTris.push_back({
                vertGrid[i][j],
                vertGrid[i + 1][j],
                vertGrid[i + 1][j + 1]
                });

            // Triangolo "destro"
            outTris.push_back({
                vertGrid[i][j],
                vertGrid[i + 1][j + 1],
                vertGrid[i][j + 1]
                });
        }
        // Triangolo in fondo alla riga
        outTris.push_back({
            vertGrid[i][i],
            vertGrid[i + 1][i],
            vertGrid[i + 1][i + 1]
            });
    }
}

void subdivideGeometry_T(
    const std::vector<Eigen::Vector3d>& baseVerts,
    const std::vector<std::array<unsigned, 3>>& baseTris,
    int b,
    std::vector<Eigen::Vector3d>& outVerts,
    std::vector<std::vector<unsigned>>& outTris )
{
    std::map<QuantizedVector, unsigned> vertexMap;
	double epsilon = 1e-6;

    for (const auto& tri : baseTris) {
        Eigen::Vector3d A = baseVerts[tri[0]];
        Eigen::Vector3d B = baseVerts[tri[1]];
        Eigen::Vector3d C = baseVerts[tri[2]];

        subdivideOnce_T(A, B, C, b, outVerts, outTris, vertexMap, epsilon );
    }
}



// funzione che computa il duale, prende in ingresso i vertici e le facce iniziali e costruisce vertici e facce del duale
void computeDualMesh(
    const std::vector<Eigen::Vector3d>& verts,
    const std::vector<std::vector<unsigned>>& tris,
    std::vector<Eigen::Vector3d>& dualVerts,
    std::vector<std::vector<unsigned>>& dualFaces)
{
    dualVerts.clear();
    dualFaces.clear();

    // valuta il centro di ogni faccia triangolare
    std::vector<Eigen::Vector3d> faceCenters;
    for (const auto& tri : tris) {
        if (tri.size() != 3) continue;
        Eigen::Vector3d center = (verts[tri[0]] + verts[tri[1]] + verts[tri[2]]) / 3.0;
        center.normalize();
        faceCenters.push_back(center);
    }

    dualVerts = faceCenters;

    // mappa i vertici originali alle facce adiacenti
    std::unordered_map<unsigned, std::vector<unsigned>> vertexToFaces;
    for (unsigned f = 0; f < tris.size(); ++f) {
		const auto& tri = tris[f];
		if (tri.size() != 3) continue; // salta facce non triangolari
        for (unsigned ve : tri) {
            vertexToFaces[ve].push_back(f);
        }
    }

    // costruisce le facce del duale
    for (unsigned v = 0; v < verts.size(); ++v) {
        const auto& faceIds = vertexToFaces[v];
        if (faceIds.size() < 3) continue;


		// calcola i vettori ortogonali e gli angoli per ordinare le facce
        Eigen::Vector3d center = verts[v];
        center.normalize();
		Eigen::Vector3d ref = (faceCenters[faceIds[0]] - center).normalized(); // usa il primo centro come riferimento
		Eigen::Vector3d orto = center.cross(ref).normalized(); // vettore ortogonale al primo centro

        
        std::vector<std::pair<double, unsigned>> angles;
        for (unsigned f : faceIds) {
            Eigen::Vector3d dir = (faceCenters[f] - center).normalized();
            double angle = std::atan2(dir.dot(orto), dir.dot(ref)); // calcola l'angolo rispetto alla direzione x (ref)
            angles.emplace_back(angle, f);
        }

        std::sort(angles.begin(), angles.end());

        std::vector<unsigned> ordered;
        for (auto& [angle, f] : angles) {
            ordered.push_back(f);
        }


        dualFaces.push_back(ordered);
    }
}




