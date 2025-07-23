# pragma once

#include<iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <map>
#include <gtest/gtest.h>
#include <GeodesicGenerator.hpp>
#include <GeometryUtils.hpp>
#include <ShortestPath.hpp>
#include <PolygonalMesh.hpp>


namespace PolygonalLibrary {

	TEST(TestGeometry, TestNumDual)   // test numerosità nel duale
	{
		unsigned q = 3;   
		array p = { 4,5 };
		array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		for (unsigned b_val : b) {
			for (unsigned p_val : p) {
				PolygonalMesh mesh = generateGeodesicMesh(p_val, q, b_val);
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
		array q = { 3,4,5 };
		array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		for (unsigned b_val : b) {
			for (unsigned q_val : q) {
				PolygonalMesh mesh = generateGeodesicMesh(p, q_val, b_val);
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

	TEST(TestGeometry, TestVerticesNormalization)   // test correttezza normalizzazione
	{
		PolygonalMesh mesh = generateGeodesicMesh(3, 5, 5);
		for (unsigned i = 0; i < mesh.NumVertices; ++i) {
			Eigen::Vector3d v = mesh.Cell0DsCoordinates.row(i);
			double vnorm = v.norm();
			EXPECT_NEAR(vnorm, 1.0, 1e-5);
		}
	}

}