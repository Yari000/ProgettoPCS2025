#pragma once

#include "PolygonalMesh.hpp"
#include "GeometryUtils.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>

namespace PolygonalLibrary {

	PolygonalMesh generateGeodesicMesh(unsigned int p, unsigned int q, unsigned int b,
	unsigned start_v=-1, unsigned end_v=-1);

}


