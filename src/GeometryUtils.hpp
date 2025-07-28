#pragma once

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <cmath>


struct QuantizedVector {
	int x, y, z;

	QuantizedVector(const Eigen::Vector3d& v, double epsilon = 1e-4) {
		x = int (std::round(v.x() / epsilon));
		y = int (std::round(v.y() / epsilon));
		z = int (std::round(v.z() / epsilon));
	}

	bool operator<(const QuantizedVector& other) const {
		if (x != other.x) return x < other.x;
		if (y != other.y) return y < other.y;
		return z < other.z;
	}
};


void generateTetrahedron(std::vector<Eigen::Vector3d>& vertices, std::vector<std::array<unsigned, 3>>& faces);
void generateOctahedron(std::vector<Eigen::Vector3d>& vertices, std::vector<std::array<unsigned, 3>>& faces);
void generateIcosahedron(std::vector<Eigen::Vector3d>& vertices, std::vector<std::array<unsigned, 3>>& faces);

void subdivideOnce_T(
	const Eigen::Vector3d& A,
	const Eigen::Vector3d& B,
	const Eigen::Vector3d& C,
	int b,
	std::vector<Eigen::Vector3d>& outVerts,
	std::vector<std::vector<unsigned>>& outTris,
	std::map<QuantizedVector, unsigned>& vertexMap,
    double epsilon=1e-6);

void subdivideGeometry_T(
	const std::vector<Eigen::Vector3d>& baseVerts,
	const std::vector<std::array<unsigned, 3>>& baseTris,
	int b,
	std::vector<Eigen::Vector3d>& outVerts,
	std::vector<std::vector<unsigned>>& outTris);

void computeDualMesh(
	const std::vector<Eigen::Vector3d>& verts,
	const std::vector<std::vector<unsigned>>& tris,
	std::vector<Eigen::Vector3d>& dualVerts,
	std::vector<std::vector<unsigned>>& dualFaces);

void rotateTriangle(
	const Eigen::Vector3d& A,
	const Eigen::Vector3d& B,
	const Eigen::Vector3d& C,
	int orientation,
	Eigen::Vector3d& Arot,
	Eigen::Vector3d& Brot,
	Eigen::Vector3d& Crot);


void subdivideOnce_T2(
	const Eigen::Vector3d& A,
	const Eigen::Vector3d& B,
	const Eigen::Vector3d& C,
	int b,
	std::vector<Eigen::Vector3d>& outVerts,
	std::vector<std::vector<unsigned>>& outTris,
	std::map<QuantizedVector, unsigned>& vertexMap,
	double epsilon = 1e-6);

void subdivideGeometry_T2(
	const std::vector<Eigen::Vector3d>& baseVerts,
	const std::vector<std::array<unsigned, 3>>& baseTris,
	int b,
	std::vector<Eigen::Vector3d>& outVerts,
	std::vector<std::vector<unsigned>>& outTris);



