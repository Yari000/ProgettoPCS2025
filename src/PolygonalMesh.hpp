#pragma once

#include <iostream>
#include "Eigen/Eigen"


namespace PolygonalLibrary {

struct PolygonalMesh
{
	unsigned int NumCell0Ds = 0;                 // numero di vertici della mesh
	std::vector<unsigned int> Cell0DsId = {};    // id dei vertici della mesh
	Eigen::MatrixXd Cell0DsCoordinates = {};     // matrice delle coordinate dei vertici della mesh
    
	unsigned int NumCell1Ds = 0;                 // numero di lati della mesh
	std::vector<unsigned int> Cell1DsId = {};    // id dei lati della mesh
	Eigen::MatrixXi Cell1DsExtrema = {};         // lista degli id dei vertici per ogni lato della mesh
    
	unsigned int NumCell2Ds = 0;                 // numero di facce della mesh
	std::vector<unsigned int> Cell2DsId = {};    // id delle facce della mesh
	std::vector<unsigned int> NumVerticesPerCell2D = {};             // numero di vertici per ogni faccia della mesh
	std::vector<unsigned int> NumEdgesPerCell2D = {};                // numero di lati per ogni faccia della mesh
	std::vector<std::vector<unsigned int>> Cell2DsVertices = {};     // lista di vertici di ogni cella 2D
	std::vector<std::vector<unsigned int>> Cell2DsEdges = {};        // lista dei lati di ogni cella 2D

	unsigned int NumCell3Ds = 0;                 // numero di celle 3D della mesh (1)
	std::vector<unsigned int> Cell3DsId = {};    // id delle celle 3D della mesh  (1)
	std::vector<std::vector<unsigned int>> Cell3DsVertices = {};     // liste dei vertici, lati e facce di ogni cella 3D    
    std::vector<std::vector<unsigned int>> Cell3DsEdges = {}; 
    std::vector<std::vector<unsigned int>> Cell3DsFaces = {};


	unsigned int NumVertices = 0;                // numero di vertici della mesh
	unsigned int NumEdges = 0;                   // numero di lati della mesh
	unsigned int NumFaces = 0;                   // numero di facce della mesh
};

}
