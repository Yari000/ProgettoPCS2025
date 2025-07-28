# pragma once

#include<iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <map>
#include <gtest/gtest.h>
#include <GeodesicGenerator.hpp>
#include <GeometryUtils.hpp>
#include <PolygonalMesh.hpp>
#include <ShortestPath.hpp>


namespace PolygonalLibrary {

	TEST(TestGeometry, TestNumDual)   // test numerosità nel duale
	{
		unsigned q = 3;   
		std::array p = { 4,5 };
		std::array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		for (unsigned b_val : b) {
			for (unsigned p_val : p) {
				PolygonalMesh mesh = generateGeodesicMesh(p_val, q, b_val,0);
				if (p_val == 4) {
					EXPECT_EQ(mesh.NumVertices, 8*b_val*b_val);
					EXPECT_EQ(mesh.NumFaces, 4*b_val*b_val +2);
					EXPECT_EQ(mesh.NumEdges, 12*b_val*b_val);
				}
				else {
					EXPECT_EQ(mesh.NumVertices, 20*b_val*b_val);
					EXPECT_EQ(mesh.NumFaces, 10*b_val*b_val +2);
					EXPECT_EQ(mesh.NumEdges, 30*b_val*b_val);
				}
				
			}
			
		}

	}


	TEST(TestGeometry, TestNum)  // test numerosità della mesh
	{
		unsigned p = 3;
		std::array q = { 3,4,5 };
		std::array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		for (unsigned b_val : b) {
			for (unsigned q_val : q) {
				PolygonalMesh mesh = generateGeodesicMesh(p, q_val, b_val, 0);
				if (q_val == 3) {
					EXPECT_EQ(mesh.NumVertices, 2*b_val*b_val +2);
					EXPECT_EQ(mesh.NumFaces, 4*b_val*b_val );
					EXPECT_EQ(mesh.NumEdges, 6*b_val*b_val);
				}
				else if (q_val == 4) {
					EXPECT_EQ(mesh.NumVertices, 4 * b_val * b_val +2);
					EXPECT_EQ(mesh.NumFaces, 8 * b_val * b_val );
					EXPECT_EQ(mesh.NumEdges, 12 * b_val * b_val);
				}
				else if (q_val == 5) {
					EXPECT_EQ(mesh.NumVertices, 10 * b_val * b_val +2);
					EXPECT_EQ(mesh.NumFaces, 20 * b_val * b_val);
					EXPECT_EQ(mesh.NumEdges, 30 * b_val * b_val);
				}
			}
		}
	}

	TEST(TestGeometry, TestNum2)
		{
			unsigned p = 3;
	        std::array q = { 3,4,5 };
			std::array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
			for (unsigned b_val : b) {
				for (unsigned q_val : q) {
					unsigned T = 3*b_val * b_val;
					PolygonalMesh mesh = generateGeodesicMesh(p, q_val, b_val, b_val);
					if (q_val == 3) {
						EXPECT_EQ(mesh.NumVertices, 2 * T + 2) << "Errore per q=3, b=" << b_val;
						EXPECT_EQ(mesh.NumFaces, 4 * T);
						EXPECT_EQ(mesh.NumEdges, 6 * T);
					}
					else if (q_val == 4) {
						EXPECT_EQ(mesh.NumVertices, 4 * T + 2) << "Errore per q=4, b=" << b_val;
						EXPECT_EQ(mesh.NumFaces, 8 * T);
						EXPECT_EQ(mesh.NumEdges, 12 * T);
					}
					else if (q_val == 5) {
						EXPECT_EQ(mesh.NumVertices, 10 * T + 2) << "Errore per q=5, b=" << b_val;
						EXPECT_EQ(mesh.NumFaces, 20 * T);
						EXPECT_EQ(mesh.NumEdges, 30 * T);
					}
				}
			}








		}

	TEST(TestGeometry, TestVerticesNormalization)   // test correttezza normalizzazione
	{
		PolygonalMesh mesh = generateGeodesicMesh(3, 5, 5,0);
		for (unsigned i = 0; i < mesh.NumVertices; ++i) {
			Eigen::Vector3d v = mesh.Cell0DsCoordinates.row(i);
			double vnorm = v.norm();
			EXPECT_NEAR(vnorm, 1.0, 1e-5);
		}
	}

	TEST(TestGeometry, EdgeConnection) {     //ogni lato ha esattamente 2 vertici distinti
		PolygonalMesh mesh = generateGeodesicMesh(3, 4, 5,0);
		for (unsigned i = 0; i < mesh.NumEdges; ++i) {
			unsigned v1 = mesh.Cell1DsExtrema(i, 0);
			unsigned v2 = mesh.Cell1DsExtrema(i, 1);
			EXPECT_LT(v1, mesh.NumVertices);
			EXPECT_LT(v2, mesh.NumVertices);
			EXPECT_NE(v1, v2);
		}
	}

	TEST(TestGeometry, ZeroPath) {    //test percorso nullo
		PolygonalMesh mesh = generateGeodesicMesh(5, 3, 2, 2);
	std::vector<bool> vertexOnPath, edgeOnPath;
		double len = 0.0;
		unsigned int numedge = 0;
		PolygonalLibrary::computeShortestPath(mesh, 2, 2, vertexOnPath, edgeOnPath, len, numedge);
		EXPECT_EQ(len, 0.0);
		EXPECT_EQ(numedge, 0);
		EXPECT_TRUE(vertexOnPath[2]);
	}

	static std::vector<unsigned> vertexVal(const PolygonalMesh& mesh) {
		std::vector<unsigned> valences(mesh.NumVertices, 0);
		for (unsigned e = 0; e < mesh.NumEdges; ++e) {
			unsigned v1 = mesh.Cell1DsExtrema(e, 0);
			unsigned v2 = mesh.Cell1DsExtrema(e, 1);
			if (v1 < mesh.NumVertices && v2 < mesh.NumVertices) {
				valences[v1]++;
				valences[v2]++;
			}
		}
		return valences;
	}

	TEST(TestGeometry, VertexValence) {
		struct Case { unsigned int q; unsigned int specCount; unsigned int specVal; };
		std::vector<Case> cases = {
			{3, 4, 3},   // tetraedro
			{4, 6, 4},   // ottaedro
			{5, 12, 5}   // icosaedro
		};

		unsigned int b = 3;
		unsigned int p = 3;

		for (const auto& c : cases) {
			PolygonalMesh mesh = generateGeodesicMesh(p, c.q, b, 0);
			auto valences = vertexVal(mesh);

			unsigned total = mesh.NumVertices;
			unsigned exspecial = c.specCount;
			unsigned exregular = total - exspecial;
			unsigned countSpecial = std::count(valences.begin(), valences.end(), c.specVal);
			unsigned countRegular = std::count(valences.begin(), valences.end(), 6u);

			EXPECT_EQ(countSpecial, exspecial) << "Errore nel numero di vertici 'speciali' per q=" << c.q;
			EXPECT_EQ(countRegular, exregular) << "Errore nel numero di vertici 'regolari' per q=" << c.q;

			for (unsigned v : valences) {
				EXPECT_TRUE(v == c.specVal || v == 6u) << "Valenza non valida: " << v << " per q=" << c.q;
			}
		}
	}

	TEST(TestGeometry, FaceConsistency) {
		PolygonalMesh mesh = generateGeodesicMesh(3, 5, 3, 0);

		std::map<std::pair<unsigned, unsigned>, unsigned int> edgefacecount;

		for (unsigned int f = 0; f < mesh.NumFaces; ++f) {
			const auto& verts = mesh.Cell2DsVertices[f];
			const auto& edges = mesh.Cell2DsEdges[f];

			EXPECT_EQ(verts.size(), edges.size()) << "Numero di vertici e lati non corrispondono nella faccia " << f;

			unsigned n = verts.size();
			for (unsigned int i = 0; i < n; ++i) {
				unsigned v1 = verts[i];
				unsigned v2 = verts[(i + 1) % n];
				bool found = false;
				for (unsigned edgeId : edges) {
					unsigned int a = mesh.Cell1DsExtrema(edgeId, 0);
					unsigned int b = mesh.Cell1DsExtrema(edgeId, 1);
					if ((a == v1 && b == v2) || (a == v2 && b == v1)) {
						found = true;
						break;
					}
				}
				EXPECT_TRUE(found) << "Lato non trovato per vertici " << v1 << " e " << v2 << " nella faccia " << f;
			}

			for (unsigned j = 0; j < n; ++j) {
				unsigned v1 = std::min(verts[j], verts[(j + 1) % n]);
				unsigned v2 = std::max(verts[j], verts[(j + 1) % n]);
				edgefacecount[{v1, v2}]++;
			}
		}

		for (auto& k : edgefacecount) {
			EXPECT_LE(k.second, 2u) << "Lato " << k.first.first << "-" << k.first.second << " appare in più di due facce (impossibile)";
		}
	}

}