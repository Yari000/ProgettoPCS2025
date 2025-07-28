#pragma once

#include "PolygonalMesh.hpp"
#include <vector>
#include <Eigen/Dense>


namespace PolygonalLibrary {
	void computeShortestPath(
		PolygonalMesh& mesh,
		unsigned int startVertex,
		unsigned int endVertex,
		std::vector<bool>& vertexOnPath,
		std::vector<bool>& edgeOnPath,
	    double& pathLength,
		unsigned& numedgeOnPath);

}

double edgeLength(const Eigen::MatrixXd& coords, int v1, int v2);
