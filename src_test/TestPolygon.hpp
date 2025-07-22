# pragma once

#include<iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <map>
#include <gtest/gtest.h>
#include <GeodesicGenerator.hpp>
#include <GeometryUtils.hpp>

namespace PolygonalLibrary {

	TEST(TestGeometry, TestNumVer)
	{
		unsigned p = 4;
		unsigned q = 3;
		array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		for (unsigned b_val : b) {
			PolygonalMesh mesh = generateGeodesicMesh(p, q, b_val);
			EXPECT_EQ(mesh.NumVertices, 8*b_val*b_val);
		}

	}


	TEST(TestGeometry, TestNumFaces)
	{
		unsigned p = 3;
		unsigned q = 5;
		array b = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		for (unsigned b_val : b) {
			PolygonalMesh mesh = generateGeodesicMesh(p, q, b_val);
			EXPECT_EQ(mesh.NumFaces, 20 * b_val * b_val);
		}
	}

}