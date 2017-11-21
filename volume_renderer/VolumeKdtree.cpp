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
	rootBox.max[0] = (float)x;
	rootBox.max[1] = (float)y;
	rootBox.max[2] = (float)z;
	rootBox.volMax[0] = rootBox.max[0] - 1;
	rootBox.volMax[1] = rootBox.max[1] - 1;
	rootBox.volMax[2] = rootBox.max[2] - 1;

	// Analyze input data dimensions to define size of kd-tree
	int numXsplits = (int)(log(X) / log(2));
	int numYsplits = (int)(log(Y) / log(2));
	int numZsplits = (int)(log(Z) / log(2));
	//oneCellDepth = numXsplits + numYsplits + numZsplits; // numTotalSplits (not counting root)
	treeDepth = numXsplits + numYsplits + numZsplits; // total number of splits (not counting root)
	
	distanceMap.resize(treeDepth + 1, 0); // treeDepth + root
	distanceSums.resize(treeDepth + 1, 0.0);
	nonzeroNodeCounts.resize(treeDepth + 1, 0.0);
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
	compressTree(true);

	// print info
	std::cout << "Tree Size: " << (double)tree.bytes() / 1e9 << " GB" << std::endl;
	for (int i = 0; i < treeDepth + 1; i++) {
		std::cout << (int)distanceMap[i] << std::endl;
	}
}

void VolumeKdtree::buildRecursive(int64_t idx, int depth, BoundingBox box, byte parentScalar) {

	tempTree[idx].parentScalar = parentScalar;

	// Loop through volume cells within bounding box & sum all scalar values
	double sumScalar = 0;
	for (int64_t x = box.volMin[0]; x < box.volMax[0] + 1; x++) {
		for (int64_t y = box.volMin[1]; y < box.volMax[1] + 1; y++) {
			for (int64_t z = box.volMin[2]; z < box.volMax[2] + 1; z++) {
				sumScalar += (double)(*data)[getCell(x, y, z)];
			}
		}
	}

	// Compute average scalar value of all volume cells
	tempTree[idx].scalar = (byte)(sumScalar / (double)box.getNumCells());

	// 1st pass of compression: record delta scalar stat
	byte estScalar;
	encodeDeltaScalar(idx, depth, parentScalar, estScalar, false);

	// Recurse on children if not at maximum tree depth
	if (depth < treeDepth) {
		int splitDim = depth % MAX_DIM;
		Point3f extent = box.getExtents();

		// try to achieve full resolution first
		if (extent.prod() > 1.0f) {
			while (extent[splitDim] == 1) {
				splitDim = splitDim == 2 ? 0 : splitDim + 1;
			}
		}
		
		float thisMid = (box.min[splitDim] + box.max[splitDim]) / 2.0f;
		float thisMax = box.max[splitDim];

		// left child
		box.max[splitDim] = thisMid;
		if (extent[splitDim] > 1)
			box.volMax[splitDim] = box.max[splitDim] - 1;
		buildRecursive((2 * idx) + 1, depth + 1, box, estScalar);

		// right child
		box.min[splitDim] = thisMid;
		box.max[splitDim] = thisMax;
		if (extent[splitDim] > 1) {
			box.volMin[splitDim] = thisMid;
			box.volMax[splitDim] = thisMax - 1;
		}
		buildRecursive((2 * idx) + 2, depth + 1, box, estScalar);
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

void VolumeKdtree::compressTree(bool computeMap) {

	// Populate distanceMap values from temporary vectors
	if (computeMap) {
		for (int depth = 0; depth <= treeDepth; depth++) {
			distanceMap[depth] = (byte)(distanceSums[depth] / nonzeroNodeCounts[depth]);
		}
	}
	
	// Compress nodes in depth first order,starting from root
	numZeroNodes = 0;
	compressTreeRecursive(0, 0);
	std::cout << "NUMBER ZERO: " << numZeroNodes << std::endl;

	// average the mse
	double numLeaves = pow(2.0, (double)treeDepth);
	mse = mse / numLeaves;
	std::cout << "MSE: " << mse << std::endl;
	std::cout << "Max Error: " << maxError << std::endl;

	//if (!added) {
	//	addLevels(3);
		//added = true;
	//}

	//pruneTree();

	// deallocate space for temporary vectors
	std::vector<TempNode>().swap(tempTree);
	std::vector<double>().swap(nonzeroNodeCounts);
	std::vector<double>().swap(distanceSums);
}


void VolumeKdtree::pruneTree() {
	// loop thru leaf nodes
	int64_t numNodes = (int64_t)(pow(2.0, (double)treeDepth + 1.0) - 1.0);
	int64_t numLeaves = pow(2, treeDepth);
	int64_t numConverted = 0;

	for (int64_t i = numNodes - numLeaves; i < numNodes; i++) {
		// get code of leaf node
		int64_t nodeIdx = i;

		// go to parent while code == 0, turning codes to 3
		while (tree[nodeIdx] == 0) {
			tree[nodeIdx] = 0; // CHANGE TO 3!
			nodeIdx = (nodeIdx - 1) / 2; // index of parent node
			numConverted++;
		}

	}

	std::cout << "NUMBER CONVERTED: " << numConverted << std::endl;
}
/*
void VolumeKdtree::addLevels(int numLevels) {

	// increase tree depth
	int oldDepth = treeDepth;
	treeDepth += numLevels;
	std::cout << "NEW DEPTH: " << treeDepth << std::endl;

	// resize vectors
	distanceMap.resize(treeDepth + 1, 0); // treeDepth + root
	distanceSums.resize(treeDepth + 1, 0.0);
	nonzeroNodeCounts.resize(treeDepth + 1, 0.0);
	int64_t numNodes = (int64_t)(pow(2.0, (double)treeDepth + 1.0) - 1.0);
	tempTree.resize(numNodes);
	tree.resize(numNodes);

	// start building from leaves of trees
	buildFromLeaves(0, 0, 0, rootBox, oldDepth);
	for (int i = 0; i < numLevels; i++) {
		distanceMap[oldDepth + i + 1] = pow(2,7-i);
	}

	mse = 0;
	maxError = 0;
	
	// Compress the temporary tree from leaves
	compressTree(false);
}
*/
void VolumeKdtree::buildFromLeaves(int64_t idx, int depth, byte parentScalar, BoundingBox box, int leafDepth) {

	// Comput scalar for this node
	byte scalar = parentScalar;
	int code = tree[idx];
	if (code == 1)
		scalar = std::min(scalar + distanceMap[depth], 255);
	else if (code == 2)
		scalar = std::max(scalar - distanceMap[depth], 0);

	// Subdivide node
	int splitDim = depth % MAX_DIM;
	Point3f extent = box.getExtents();

	if (extent.prod() > 1.0f) {
		while (extent[splitDim] == 1) {
			splitDim = splitDim == 2 ? 0 : splitDim + 1;
		}
	}

	float thisMid = (box.min[splitDim] + box.max[splitDim]) / 2.0f;
	float thisMax = box.max[splitDim];

	// left child
	box.max[splitDim] = thisMid;
	if (extent[splitDim] > 1)
		box.volMax[splitDim] = box.max[splitDim] - 1;

	if (depth == leafDepth) {
		buildRecursive(2 * idx + 1, depth + 1, box, scalar);
	}
	else {
		buildFromLeaves(2 * idx + 1, depth + 1, scalar, box, leafDepth);
	}

	// right child
	box.min[splitDim] = thisMid;
	box.max[splitDim] = thisMax;
	if (extent[splitDim] > 1) {
		box.volMin[splitDim] = thisMid;
		box.volMax[splitDim] = thisMax - 1;
	}
	
	if (depth == leafDepth) {
		buildRecursive(2 * idx + 2, depth + 1, box, scalar);
	}
	else {
		buildFromLeaves(2 * idx + 2, depth + 1, scalar, box, leafDepth);
	}
}

void VolumeKdtree::compressTreeRecursive(int64_t idx, int depth) {

	

	// Find scalar value & code for current node
	byte scalar;
	byte code = encodeDeltaScalar(idx, depth, tempTree[idx].parentScalar, scalar, true);
	tree[idx] = (int)code;
	if (tree[idx] == 0)
		numZeroNodes++;

	//if (depth > 23 && tempTree[idx].parentScalar - tempTree[idx].scalar > 32) {
	//	added = true;
	//	std::cout << (int)tempTree[idx].scalar << " " << (int)tempTree[idx].parentScalar << " " << (int)code << " " << (int)distanceMap[depth] << std::endl;
	//}

	// Measure mse of leaf nodes
	if (depth == treeDepth) {
		double error = abs((double)scalar - (double)tempTree[idx].scalar);
		mse += pow(error, 2);
		if (error > maxError)
			maxError = error;
	}

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
	int64_t numOutputNodes = pow(2, cutDepth);
	outData.resize(numOutputNodes);
	std::vector<double> sumData(numOutputNodes);
	std::vector<double> countData(numOutputNodes);
	levelCutRecursive(0, 0, 0, rootBox, cutDepth, sumData, countData);
	std::transform(sumData.begin(), sumData.end(), countData.begin(), outData.begin(), std::divides<double>());
}

void VolumeKdtree::levelCutRecursive(int64_t idx, int depth, byte parentScalar,
	BoundingBox box, int cutDepth, std::vector<double> &sumData, std::vector<double> &countData) {

	// For easy reference
	int code = tree[idx];
	
	// Compute scalar value for this node
	byte scalar = parentScalar;
	if (code == 1)
		scalar = std::min(scalar + distanceMap[depth], 255);
	else if (code == 2)
		scalar = std::max(scalar - distanceMap[depth], 0);

	// Populate output vectors if at cut depth or end of tree branch
	//if (code == 3 || depth == cutDepth) {
	if (depth == cutDepth) {
		//box.print();
		for (int64_t x = box.volMin[0]; x <= box.volMax[0]; x++) {
			for (int64_t y = box.volMin[1]; y <= box.volMax[1]; y++) {
				for (int64_t z = box.volMin[2]; z <= box.volMax[2]; z++) {
					int64_t c = getCell(x, y, z);
					//std::cout << c << std::endl;
					sumData[c] += scalar;
					countData[c]++;
				}
			}
		}
		return;
	}

	// OTHERWISE, recurse
	int splitDim = depth % MAX_DIM;
	Point3f extent = box.getExtents();

	if (extent.prod() > 1.0f) {
		while (extent[splitDim] == 1) {
			splitDim = splitDim == 2 ? 0 : splitDim + 1;
		}
	}

	float thisMid = (box.min[splitDim] + box.max[splitDim]) / 2.0f;
	float thisMax = box.max[splitDim];

	// left child
	box.max[splitDim] = thisMid;
	if (extent[splitDim] > 1)
		box.volMax[splitDim] = box.max[splitDim] - 1;
	levelCutRecursive(2 * idx + 1, depth + 1, scalar, box, cutDepth, sumData, countData);

	// right child
	box.min[splitDim] = thisMid;
	box.max[splitDim] = thisMax;
	if (extent[splitDim] > 1) {
		box.volMin[splitDim] = thisMid;
		box.volMax[splitDim] = thisMax - 1;
	}
	levelCutRecursive(2 * idx + 2, depth + 1, scalar, box, cutDepth, sumData, countData);
}

void VolumeKdtree::save(std::string filename) {
	// Get size of  KDOctree 
	int64_t treeSize = tree.bytes();

	// Write out the data
	std::ofstream out(filename, std::ios::out | std::ios::binary);

	out.write(reinterpret_cast<char *>(&rootBox), sizeof(BoundingBox));
	out.write(reinterpret_cast<char *>(&treeDepth), sizeof(int));
	out.write(reinterpret_cast<char *>(&X), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Z), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&distanceMap[0]), treeDepth + 1);
	out.write(reinterpret_cast<char *>(&tree.bits[0]), treeSize);
	
	out.close();

	// Print info
	int64_t dataSize = treeSize + sizeof(BoundingBox) + sizeof(int) + treeDepth + 1 + (3*sizeof(int64_t));
	std::cout << "\nNew file saved: " << filename << std::endl;
	std::cout << "File size: " << (double)dataSize / 1e9 << " GB" << std::endl;
}

void VolumeKdtree::open(std::string filename) {

	// Open the stream
	std::ifstream is(filename, std::ios::in | std::ios::binary);

	// Check that the file exists
	if (!is.good()) {
		std::cout << "\nERROR! Tree File does not exist: " << filename << std::endl;
		std::cout << "Enter to exit." << std::endl;
		std::cin.ignore();
		exit(-1);
	}

	// Get tree size from file size
	is.seekg(0, is.end);
	int64_t fileSize = is.tellg();
	is.seekg(0, is.beg);

	is.read(reinterpret_cast<char *>(&rootBox), sizeof(BoundingBox));
	is.read(reinterpret_cast<char *>(&treeDepth), sizeof(int));
	is.read(reinterpret_cast<char *>(&X), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&Z), sizeof(int64_t));

	int64_t treeSize = fileSize - (sizeof(BoundingBox) + sizeof(int) + treeDepth + 1 + (3 * sizeof(int64_t)));
	std::cout << treeSize << std::endl;

	// Allocate memory for distanceMap & tree
	distanceMap.resize(treeDepth + 1);
	tree.bits.resize(treeSize);

	is.read(reinterpret_cast<char *>(&distanceMap[0]), treeDepth + 1);
	is.read(reinterpret_cast<char *>(&tree.bits[0]), treeSize);

	// Close the file
	is.close();

}