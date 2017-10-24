#include "VolumeKdtree.h"


inline int64_t VolumeKdtree::getCell(int64_t x, int64_t y, int64_t z) {
	return x + (X*y) + (X*Y*z);
}


void VolumeKdtree::build(std::vector<byte> &inData, int64_t x, int64_t y, int64_t z) {

	// Assign values to data-dependent class members
	
	data = &inData;
	X = x;
	Y = y;
	Z = z;
	rootBox.max[0] = x;
	rootBox.max[1] = y;
	rootBox.max[2] = z;

	// Analyze input data dimensions to define size of kd-tree
	int numXsplits = (int)(log(X) / log(2));
	int numYsplits = (int)(log(Y) / log(2));
	int numZsplits = (int)(log(Z) / log(2));
	treeDepth = numXsplits + numYsplits + numZsplits;
	
	
	distanceMap.resize(treeDepth + 1,0);
	distanceSums.resize(treeDepth + 1,0.0);
	nonzeroNodeCounts.resize(treeDepth + 1,0.0);
	int64_t numNodes = (int64_t)(pow(2.0, (double)treeDepth + 1.0) - 1.0);
	tempTree.resize(numNodes);
	tree.resize(numNodes);

	// Print tree info
	std::cout << "Tree depth: " << treeDepth << std::endl;
	std::cout << "Number of nodes: " << numNodes << std::endl;
	std::cout << "Estimated Tree Size : " << (numNodes * (2.0 / 8.0)) / 1e9 << " GB" << std::endl;

	// Call recursive build function from root - builds tempTree
	buildRecursive((int64_t)0, 0, rootBox, 0);

	// Compress the temporary tree
	compressTree();

	std::cout << "Tree Size: " << (double)tree.size() / 1e9 << " GB" << std::endl;

	for (int i = 0; i < treeDepth; i++) {
		std::cout << (int)distanceMap[i] << std::endl;
	}
}

void VolumeKdtree::buildRecursive(int64_t idx, int depth, BoundingBox box, byte parentScalar) {

	// Fill in known fields
	tempTree[idx].parentScalar = parentScalar;

	// Compute average of scalar value in node by looping thru cells in bounding box
	int64_t numCells = box.getVolume();
	double sumScalar = 0;
	int64_t cellIdx;
	for (int64_t x = box.min[0]; x < box.max[0]; x++) {
		for (int64_t y = box.min[1]; y < box.max[1]; y++) {
			for (int64_t z = box.min[2]; z < box.max[2]; z++) {
				cellIdx = getCell(x, y, z);
				sumScalar += (double)(*data)[cellIdx];
			}
		}
	}
	tempTree[idx].scalar = (byte)(sumScalar / numCells);;

	//if (idx <= 7)
	//	std::cout << idx << " " << (int)tempTree[idx].parentScalar << " " << (int)tempTree[idx].scalar << std::endl;

	// 1st pass of compression: record delta scalar stat
	byte estScalar;
	encodeDeltaScalar(idx, depth, parentScalar, estScalar, false);

	// Split bounding box & recurse on children
	if (depth < treeDepth) {
		int splitDim = depth % MAX_DIM;
		int3D extent;
		box.getExtent(extent);
		while (extent[splitDim] == 1) {
			splitDim = splitDim == 2 ? 0 : splitDim + 1;
		}
		std::vector<BoundingBox> children = box.split2(splitDim);
		buildRecursive((2 * idx) + 1, depth + 1, children[0], tempTree[idx].scalar); // left 
		buildRecursive((2 * idx) + 2, depth + 1, children[1], tempTree[idx].scalar); // right
	}
	

}

byte VolumeKdtree::encodeDeltaScalar(int64_t idx, int depth,
	byte parent, byte &estimate, bool useMap) {

	// For easy reference
	byte truth = tempTree[idx].scalar;

	// Compute sclar difference between parent node and current node
	int delta = (int)truth - (int)parent;
	byte distance = (byte)abs(delta);

	// Return 0 if no change from parent
	if (delta == 0)
		return 0;

	// Define mapped delta based on deltaScalarMap or estimate from current available node data
	byte mappedDistance = useMap ? distanceMap[depth] : byte(double(distanceSums[depth] + distance) 
		/ double(nonzeroNodeCounts[depth] + 1));

	// Return the delta code with the minimum error
	byte minEstError = delta;
	byte code = 0;

	// Calculate error for each code case
	byte nochangeError = distance;
	byte plusEst = std::min(255, (int)parent + (int)mappedDistance);
	byte plusError = abs(truth - plusEst);
	byte minusEst = std::max(0, (int)parent - (int)mappedDistance);
	byte minusError = abs(truth - minusEst);

	// Return code 0 if that produces the least error
	if ((nochangeError < plusError) && (nochangeError < minusError)) {
		estimate = parent;
		return 0;
	}

	// If not returning 0, update distanceSums & nonzeroNodeCounts,
	// if distanceMap not yet populated
	if (!useMap) {
		distanceSums[depth] += (double)distance;
		nonzeroNodeCounts[depth]++;
	}

	// Return code 1 if that produces less error than code 2
	if (plusError < minusError) {
		estimate = plusEst;
		return 1;
	}

	// At this point, code 2 is the only case remaining
	estimate = minusEst;
	return 2;
}

void VolumeKdtree::compressTree() {

	// Populate distanceMap values from temporary vectors
	for (int depth = 0; depth < treeDepth; depth++) {
		distanceMap[depth] = (byte)(distanceSums[depth] / nonzeroNodeCounts[depth]);
	}

	//std::cout << (int)distanceMap[0] << " " << (int)distanceMap[1] << " " << (int)distanceMap[2] << std::endl;

	// Compress nodes in depth first order,starting from root
	compressTreeRecursive(0, 0);

	// deallocate space for temporary vectors
	std::vector<TempNode>().swap(tempTree);
	std::vector<double>().swap(nonzeroNodeCounts);
	std::vector<double>().swap(distanceSums);
}

void VolumeKdtree::compressTreeRecursive(int64_t idx, int depth) {

	// Find scalar value & code for current node
	byte scalar;
	byte code = encodeDeltaScalar(idx, depth, tempTree[idx].parentScalar, scalar, true);
	tree[idx] = (int)code;

	//if (idx <= 7)
	//	std::cout << idx << " " << (int)tree[idx] << " " << (int)scalar << std::endl;

	if (depth < treeDepth) {
		int64_t leftIdx = (2 * idx) + 1;
		int64_t rightIdx = leftIdx + 1;

		// Update parent color of node's children
		tempTree[leftIdx].parentScalar = scalar;
		tempTree[rightIdx].parentScalar = scalar;

		// Recurse on children
		compressTreeRecursive(leftIdx, depth + 1); 
		compressTreeRecursive(rightIdx, depth + 1); 
	}
}

void VolumeKdtree::levelCut(int cutDepth, std::vector<byte> &outData) {
	outData.resize((int64_t)(pow(2.0, cutDepth)));
	levelCutRecursive(0, 0, 0, rootBox, cutDepth, outData);
}

void VolumeKdtree::levelCutRecursive(int idx, int depth, byte parentScalar, 
	BoundingBox box, int cutDepth, std::vector<byte> &outData) {

	//if (depth == treeDepth)
	//	box.print();
	
	byte scalar = parentScalar;

	int code = tree[idx];

	if (code == 1)
		scalar += distanceMap[depth];
	else if (code == 2)
		scalar -= distanceMap[depth];

	if (depth < cutDepth) {
		int splitDim = depth % MAX_DIM;
		int3D extent;
		box.getExtent(extent);
		while (extent[splitDim] == 1) {
			splitDim = splitDim == 2 ? 0 : splitDim + 1;
		}
		std::vector<BoundingBox> children = box.split2(splitDim);
		levelCutRecursive(2 * idx + 1, depth + 1, scalar, children[0], cutDepth, outData);
		levelCutRecursive(2 * idx + 2, depth + 1, scalar, children[1], cutDepth, outData);
	}
	else {
		if (cutDepth == treeDepth) {
			//int3D boxCenter;
			//box.getCenter(boxCenter);
			int64_t c = getCell(box.min[0], box.min[1], box.min[2]);
			outData[c] = scalar;
		}
		//TODO: new getCell function for not whole volume
	}

}