#include "Octree.h"


inline int64_t Octree::getCell(int64_t x, int64_t y, int64_t z) {
	return x + (X*y) + (X*Y*z);
}

void Octree::build(std::vector<byte> * inData, int64_t x, int64_t y, int64_t z) {

	// Assign values to data-dependent class members
	data = inData;
	X = x;
	Y = y;
	Z = z;
	rootBox.max[0] = (float)x;
	rootBox.max[1] = (float)y;
	rootBox.max[2] = (float)z;

	// Analyze input data dimensions to define size of octree 
	
	// If input volume is perfect cube, use octree formulas to find exact number of tree nodes
	// and change tempTree container accordingly
	if ((X == Y) && (Y == Z)) {
		octreeStartDepth = 0;
		int64_t numCells = data->size();
		int depth = calcOctreeDepth(numCells);
		int64_t numNodes = calcNumOctreeNodes(depth);
		tempTree.resize(numNodes);
		std::cout << "Number of tree nodes: " << numNodes << std::endl;
	}
	else {
		// find number of nodes in "mini cubes"
		// mini cubes = perfect cubes that make up the full volume
		int64_t minDim = std::min(X, Y);
		minDim = std::min(minDim, Z);
		if (X != minDim)
			splitDims.push_back(0);
		if (Y != minDim)
			splitDims.push_back(1);
		if (Z != minDim)
			splitDims.push_back(2);

		int64_t miniCubeCells = minDim * minDim * minDim;
		int miniCubeDepth = calcOctreeDepth(miniCubeCells);
		int64_t miniCubeNodes = calcNumOctreeNodes(miniCubeDepth);

		// find number of binary tree levels needed to get to minicubes 
		int binTreeLevels = 0;
		int64_t volumeDims[] = { X, Y, Z };
		for (int i = 0; i < 3; i++) {
			binTreeLevels += (int)(log(volumeDims[i] / minDim) / log(2));
		}
		octreeStartDepth = binTreeLevels;

		// total # of nodes = # mini cubes * # nodes in mini cubes  + # binary tree levels + # of missing nodes
		// missing nodes = idx of first node in octreeStartDepth+1 - idx of last node in octreeStartDepth
		int64_t missingNodes = (((pow(2, octreeStartDepth) - 1) * 8) + 1) - (pow(2, octreeStartDepth + 1) - 1);
		int numMiniCubes = (int) pow(2, binTreeLevels);
		int64_t numBinTreeNodes = (int64_t) (pow(2, binTreeLevels + 1) - 1);
		int64_t numNodes = (numMiniCubes*miniCubeNodes) + numBinTreeNodes + missingNodes;
		tempTree.resize(numNodes);
		std::cout << "Number of tree nodes: " << numNodes << std::endl;
		std::cout << "Number of missing nodes: " << missingNodes << std::endl;
	}

	// Print tree info
	std::cout << "Estimated Tree Size : " << (tempTree.size() * (2.0/8.0)) / 1e9 << " GB" << std::endl;

	// Call recursive build function from root
	buildRecursive((int64_t)0, 0, rootBox, 0);
}

void Octree::buildRecursive(int64_t idx, int depth, BoundingBox box, byte parentScalar) {

	// Fill in known fields
	tempTree[idx].parentScalar = parentScalar;
	
	// Compute average of scalar value in node by looping thru cells in bounding box
	int64_t numCells = box.getVolume();
	double sumScalar = 0;
	int64_t cellIdx;
	for (int x = box.min[0]; x < box.max[0]; x++){
		for (int y = box.min[1]; y < box.max[1]; y++){
			for (int z = box.min[2]; z < box.max[2]; z++) {
				cellIdx = getCell(x, y, z);
				sumScalar += (double)(*data)[cellIdx];
			}
		}
	}
	byte avgScalar = (byte)(sumScalar / numCells);
	tempTree[idx].scalar = avgScalar;
	
	// Record delta scalar stat
	//byte estScalar;
	//encodeDeltaScalar(idx, depth, parentScalar, estScalar, false);

	// Split bounding box & recurse on children
	if (depth < octreeStartDepth) {
		// In binary part of tree, split into 2 children
		int splitDimIdx = depth % splitDims.size();
		std::vector<BoundingBox> children = box.split2(splitDims[splitDimIdx]);

		buildRecursive( (2 * idx) + 1, depth + 1, children[0], avgScalar); // left 
		buildRecursive( (2 * idx) + 2, depth + 1, children[1], avgScalar); // right
	}
	else {
		// In octree part of tree, split into 8 children
		std::vector<BoundingBox> children = box.split8();
		for (int i = 0; i < 8; i++)
			buildRecursive((8 * idx) + i + 1, depth + 1, children[i], avgScalar);
	}
}

int Octree::calcOctreeDepth(int64_t numCells) {
	return (int)(log(numCells) / log(8));
}

int64_t Octree::calcNumOctreeNodes(int depth) {
	return (int64_t)((pow(8.0, (double)depth + 1.0) - 1.0) / 7.0);
}