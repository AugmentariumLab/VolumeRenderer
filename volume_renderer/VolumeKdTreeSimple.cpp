#include "VolumeKdtreeSimple.h"


inline int64_t VolumeKdtreeSimple::getCell(int64_t x, int64_t y, int64_t z) {
	return x + (X*y) + (X*Y*z);
}

void VolumeKdtreeSimple::build() {

	// Analyze input data dimensions to define size of kd-tree
	int numXsplits = (int)(log(X) / log(2));
	int numYsplits = (int)(log(Y) / log(2));
	int numZsplits = (int)(log(Z) / log(2));
	treeDepth = numXsplits + numYsplits + numZsplits; // total number of splits (not counting root)
	distanceMap.resize(treeDepth + 1, 0); // treeDepth + root
	distanceSums.resize(treeDepth + 1, 0.0); // treeDepth + root
	distanceCounts.resize(treeDepth + 1, 0.0); // treeDepth + root
	int64_t numNodes = (int64_t)pow(2, treeDepth + 1) - 1;
	temp.resize(numNodes);
	tree.resize(numNodes);

	// Print tree info
	std::cout << "Tree depth: " << treeDepth << std::endl;
	std::cout << "Number of nodes: " << numNodes << std::endl;
	std::cout << "Estimated Tree Size : " << (numNodes * (2.0 / 8.0)) / 1e9 << " GB" << std::endl;

	// Call recursive build function starting from root
	buildRecursive(0, 0, rootMin, rootMax, 0);

	// Compute Distance Map
	for (int depth = 0; depth < treeDepth + 1; depth++) {
		distanceMap[depth] = (byte)(distanceSums[depth] / distanceCounts[depth]);
	}
	std::vector<double>().swap(distanceSums);
	std::vector<double>().swap(distanceCounts);
	

	//Fix Distance Map
	//int cycle[] = { 1,2,4,8,16,32,64,128 };
	/*
	int cycle[] = { 128,64,32,16,8,4,2,1 };
	for (int depth = 0; depth < treeDepth + 1; depth++) {
		//distanceMap[depth] = pow(2, cycle[depth % 7]);
		distanceMap[depth] = cycle[depth % 8];
	}
	*/
	
	/*
	int cycle[] = { 128,64,32,16,8,4,2,1 };
	for (int depth = treeDepth - 7; depth < treeDepth + 1; depth++) {
		//distanceMap[depth] = pow(2, cycle[depth % 7]);
		distanceMap[depth] = cycle[depth % 8];
	}
	*/

	

	// 2nd pass of compression
	compressTreeRecursive(0, 0, 0);

	// add additional levels
	for (int depth = 0; depth < treeDepth + 1; depth++) {
		std::cout << (int)distanceMap[depth] << std::endl;
	}
	//addLevels(2);
	//for (int depth = 0; depth < treeDepth + 1; depth++) {
	//	std::cout << (int)distanceMap[depth] << std::endl;
	//}


	// delete temporary tree
	std::vector<byte>().swap(temp);
}

void VolumeKdtreeSimple::buildRecursive(int64_t idx, int depth, Point3i minBound, 
	Point3i maxBound, byte parentEstimate) {

	// Loop through volume cells within bounding box & sum all scalar values
	double sumScalar = 0.0;
	double maxScalar = 0.0;
	double minScalar = 255.0;
	for (int64_t x = minBound[0]; x < maxBound[0]; x++) {
		for (int64_t y = minBound[1]; y < maxBound[1]; y++) {
			for (int64_t z = minBound[2]; z < maxBound[2]; z++) {
				sumScalar += (double)(*data)[getCell(x, y, z)];
				maxScalar = std::max(maxScalar, (double)(*data)[getCell(x, y, z)]);
				minScalar = std::min(minScalar, (double)(*data)[getCell(x, y, z)]);
			}
		}
	}

	// Compute average scalar value of all volume cells
	Point3i extent = maxBound - minBound;
	int64_t numCells = extent.prod();
	//temp[idx] = (byte)(sumScalar / numCells);
	temp[idx] = (byte)((maxScalar + minScalar) / 2.0);

	// 1st pass of compression: populate distance map
	byte scalarEstimate = encodeNode(idx, depth, parentEstimate, false);

	// Recurse on children if not at maximum tree depth
	if (depth < treeDepth) {

		// Do NOT split bounding box if only 1 cell is enclosed
		if (numCells == 1) {
			buildRecursive(2 * idx + 1, depth + 1, minBound, maxBound, scalarEstimate); // left
			buildRecursive(2 * idx + 2, depth + 1, minBound, maxBound, scalarEstimate); // right
		}

		else {
			int splitDim = depth % MAX_DIM;

			// Choose split dimension to acheive full resolution
			int i = 0;
			while (numCells > 1 && extent[splitDim] == 1) {
				splitDim = (depth + ++i) % MAX_DIM;
			}

			// Split Bounding box 
			int64_t thisMid = (minBound[splitDim] + maxBound[splitDim]) / 2;
			int64_t thisMax = maxBound[splitDim];

			// left child
			maxBound[splitDim] = thisMid;
			buildRecursive(2 * idx + 1, depth + 1, minBound, maxBound, scalarEstimate);

			// right child
			minBound[splitDim] = thisMid;
			maxBound[splitDim] = thisMax;
			buildRecursive(2 * idx + 2, depth + 1, minBound, maxBound, scalarEstimate);

		}	
	}
}

void VolumeKdtreeSimple::levelCut(int cutDepth, std::vector<byte> &outData) {

	output = &outData;
	queryDepth = cutDepth;

	// Resize output vector
	output->resize(X*Y*Z);

	// Call recursive level cut function from root
	levelCutRecursive(0, 0, rootMin, rootMax, 0);
}

int VolumeKdtreeSimple::measureMaxError() {
	int maxError = 0;
	for (int64_t i = 0; i < X*Y*Z; i++) {
		maxError = std::max((int)abs((double)(*output)[i] - (double)(*data)[i]), maxError);
	}
	return maxError;
}

double VolumeKdtreeSimple::measureMeanError() {
	double sumError = 0.0;
	int64_t count = X*Y*Z;
	for (int64_t i = 0; i < count; i++) {
		sumError += abs((double)(*output)[i] - (double)(*data)[i]);
	}
	return sumError / (double)count;
}

void VolumeKdtreeSimple::queryError(std::vector<byte> &outData) {
	output2 = &outData;
	int64_t count = X*Y*Z;
	output2->resize(count);
	for (int64_t i = 0; i < count; i++) {
		(*output2)[i] = (byte)abs((double)(*output)[i] - (double)(*data)[i]);
	}
}

void VolumeKdtreeSimple::levelCutRecursive(int64_t idx, int depth, Point3i minBound, Point3i maxBound, byte compressedParent) {

	double parentEstimate = (double)compressedParent;
	int code = (int)tree[idx];

	// Compute scalar value for this node
	byte scalar;
	if (code == 0)
		scalar = compressedParent;
	else if (code == 1)
		scalar = (byte)std::min(255.0, parentEstimate + (double)distanceMap[depth]);
	else if (code == 2)
		scalar = (byte)std::max(0.0, parentEstimate - (double)distanceMap[depth]);


	// Populate output vectors if at cut depth or end of tree branch
	if (depth == queryDepth) {
		for (int64_t x = minBound[0]; x < maxBound[0]; x++) {
			for (int64_t y = minBound[1]; y < maxBound[1]; y++) {
				for (int64_t z = minBound[2]; z < maxBound[2]; z++) {
					int64_t c = getCell(x, y, z);
					(*output)[c] = scalar;
				}
			}
		}
		return;
	}

	// Recurse on children if not at maximum tree depth
	if (depth < treeDepth) {

		Point3i extent = maxBound - minBound;
		int64_t numCells = extent.prod();

		// Do NOT split bounding box if only 1 cell is enclosed
		if (numCells == 1) {
			levelCutRecursive(2 * idx + 1, depth + 1, minBound, maxBound, scalar); // left
			levelCutRecursive(2 * idx + 2, depth + 1, minBound, maxBound, scalar); // right
		}

		else {
			int splitDim = depth % MAX_DIM;

			// Choose split dimension to acheive full resolution
			int i = 0;
			while (numCells > 1 && extent[splitDim] == 1) {
				splitDim = (depth + ++i) % MAX_DIM;
			}

			// Split Bounding box 
			int64_t thisMid = (minBound[splitDim] + maxBound[splitDim]) / 2;
			int64_t thisMax = maxBound[splitDim];

			// left child
			maxBound[splitDim] = thisMid;
			levelCutRecursive(2 * idx + 1, depth + 1, minBound, maxBound, scalar);

			// right child
			minBound[splitDim] = thisMid;
			maxBound[splitDim] = thisMax;
			levelCutRecursive(2 * idx + 2, depth + 1, minBound, maxBound, scalar);

		}
	}
}

byte VolumeKdtreeSimple::encodeNode(int64_t idx, int depth, byte compressedParent, bool useMap) {

	// For easy reference
	double parentEstimate = (double)compressedParent;
	double nodeTruth = (double)temp[idx];
	//double nodeTruth = ((double)tempMax[idx] - (double)tempMin[idx]) / 2.0;

	// Define parent distance: distance from this node's true scalar value to its compressed parent
	double parentDistance = abs(parentEstimate - nodeTruth);

	// Define master distance: distance value to be added or subtracted to this node
	double masterDistance = useMap ? (double)distanceMap[depth] : 
		(distanceSums[depth] + parentDistance) / (distanceCounts[depth] + 1.0);

	// Calculate error of *doing nothing* 
	double noneEstimate = parentEstimate;
	double noneError = parentDistance;

	// Calculate error of *adding* the master distance
	double addEstimate = std::min(255.0, parentEstimate + masterDistance);
	double addError = abs(addEstimate - nodeTruth);

	// Calculate error of *subtracting* the master difference
	double subEstimate = std::max(0.0, parentEstimate - masterDistance);
	double subError = abs(subEstimate - nodeTruth);

	// Calculate minimum error scenario
	double minError = std::min(subError, std::min(noneError, addError));

	// If doing nothing produces minimum error, set node to 0 & return corresponding estimate
	if (minError == noneError) {
		tree[idx] = 0;
		return (byte)noneEstimate;
	}

	// If adding produces minimum error, set node to 1 & return corresponding estimate
	if (minError == addError) {
		if (!useMap) {
			distanceSums[depth] += addError;
			distanceCounts[depth]++;
		}
		tree[idx] = 1;
		return (byte)addEstimate;
	}

	// If subtracting produces minimum error, set node to 2 & return corresponding estimate
	if (minError == subError) {
		if (!useMap) {
			distanceSums[depth] += subError;
			distanceCounts[depth]++;
		}
		tree[idx] = 2;
		return (byte)subEstimate;


	}
}

void VolumeKdtreeSimple::compressTreeRecursive(int64_t idx, int depth, byte compressedParent) {

	// Find scalar value & code for current node
	byte scalar = encodeNode(idx, depth, compressedParent, true);

	// Recurse on children
	if (depth < treeDepth) {
		compressTreeRecursive(2*idx + 1, depth + 1, scalar);
		compressTreeRecursive(2*idx + 2, depth + 1, scalar);
	}
}

void VolumeKdtreeSimple::save(std::string filename) {
	// Get size of  KDOctree 
	int64_t treeSize = tree.bytes();

	// Write out the data
	std::ofstream out(filename, std::ios::out | std::ios::binary);

	out.write(reinterpret_cast<char *>(&rootMin), sizeof(Point3i));
	out.write(reinterpret_cast<char *>(&rootMax), sizeof(Point3i));
	out.write(reinterpret_cast<char *>(&treeDepth), sizeof(int));
	out.write(reinterpret_cast<char *>(&X), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Z), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&distanceMap[0]), treeDepth + 1);
	out.write(reinterpret_cast<char *>(&tree.bits[0]), treeSize);

	out.close();

	// Print info
	int64_t dataSize = treeSize + (2 * sizeof(Point3i)) + sizeof(int) + treeDepth + 1 + (3 * sizeof(int64_t));
	std::cout << "\nNew file saved: " << filename << std::endl;
	std::cout << "File size: " << (double)dataSize / 1e9 << " GB" << std::endl;
}

void VolumeKdtreeSimple::open(std::string filename) {

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

	is.read(reinterpret_cast<char *>(&rootMin), sizeof(Point3i));
	is.read(reinterpret_cast<char *>(&rootMax), sizeof(Point3i));
	is.read(reinterpret_cast<char *>(&treeDepth), sizeof(int));
	is.read(reinterpret_cast<char *>(&X), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&Z), sizeof(int64_t));

	int64_t treeSize = fileSize - ((2*sizeof(Point3i)) + sizeof(int) + treeDepth + 1 + (3 * sizeof(int64_t)));
	std::cout << treeSize << std::endl;

	// Allocate memory for distanceMap & tree
	distanceMap.resize(treeDepth + 1);
	tree.bits.resize(treeSize);

	is.read(reinterpret_cast<char *>(&distanceMap[0]), treeDepth + 1);
	is.read(reinterpret_cast<char *>(&tree.bits[0]), treeSize);

	// Close the file
	is.close();

}

void VolumeKdtreeSimple::addLevels(int numLevels) {

	// increase tree depth
	int oldDepth = treeDepth;
	treeDepth += numLevels;
	std::cout << "NEW DEPTH: " << treeDepth << std::endl;

	// resize vectors
	distanceMap.resize(treeDepth + 1, 0); // treeDepth + root
	distanceSums.resize(treeDepth + 1, 0.0);
	distanceCounts.resize(treeDepth + 1, 0.0);
	int64_t numNodes = (int64_t)(pow(2.0, (double)treeDepth + 1.0) - 1.0);
	temp.resize(numNodes);
	tree.resize(numNodes);

	// start building from leaves of trees
	buildFromLeaves(0, 0, rootMin, rootMax, 0, oldDepth);

	// Compute Distance Map
	//for (int depth = 0; depth < treeDepth + 1; depth++) {
	//	distanceMap[depth] = (byte)(distanceSums[depth] / distanceCounts[depth]);
	//}
	for (int i = 0; i < numLevels; i++) {
		distanceMap[oldDepth + i + 1] = (byte)pow(2, 5 - i);
	}


	// Compress the temporary tree from leaves
	compressFromLeaves(0, 0, 0, oldDepth);
}

void VolumeKdtreeSimple::buildFromLeaves(int64_t idx, int depth, 
	Point3i minBound, Point3i maxBound, byte compressedParent, int leafDepth) {

	double parentEstimate = (double)compressedParent;
	int code = (int)tree[idx];

	// Compute scalar value for this node
	byte scalar;
	if (code == 0)
		scalar = compressedParent;
	else if (code == 1)
		scalar = (byte)std::min(255.0, parentEstimate + (double)distanceMap[depth]);
	else if (code == 2)
		scalar = (byte)std::max(0.0, parentEstimate - (double)distanceMap[depth]);

	// Recurse on children if not at maximum tree depth
	Point3i extent = maxBound - minBound;
	int64_t numCells = extent.prod();
	if (depth < treeDepth) {

		// Do NOT split bounding box if only 1 cell is enclosed
		if (numCells == 1) {
			if (depth == leafDepth) {
				buildRecursive(2 * idx + 1, depth + 1, minBound, maxBound, scalar); // left
				buildRecursive(2 * idx + 2, depth + 1, minBound, maxBound, scalar); // right
			}
			else {
				buildFromLeaves(2 * idx + 1, depth + 1, minBound, maxBound, scalar, leafDepth);
				buildFromLeaves(2 * idx + 2, depth + 1, minBound, maxBound, scalar, leafDepth);
			}
		}

		else {
			int splitDim = depth % MAX_DIM;

			// Choose split dimension to acheive full resolution
			int i = 0;
			while (numCells > 1 && extent[splitDim] == 1) {
				splitDim = (depth + ++i) % MAX_DIM;
			}

			// Split Bounding box 
			int64_t thisMid = (minBound[splitDim] + maxBound[splitDim]) / 2;
			int64_t thisMax = maxBound[splitDim];

			// left child
			maxBound[splitDim] = thisMid;
			if (depth == leafDepth)
				buildRecursive(2 * idx + 1, depth + 1, minBound, maxBound, scalar);
			else
				buildFromLeaves(2 * idx + 1, depth + 1, minBound, maxBound, scalar, leafDepth);

			// right child
			minBound[splitDim] = thisMid;
			maxBound[splitDim] = thisMax;
			if (depth == leafDepth)
				buildRecursive(2 * idx + 2, depth + 1, minBound, maxBound, scalar);
			else
				buildFromLeaves(2 * idx + 2, depth + 1, minBound, maxBound, scalar, leafDepth);

		}
	}
}

void VolumeKdtreeSimple::compressFromLeaves(int64_t idx, int depth, byte compressedParent, int leafDepth) {
	
	double parentEstimate = (double)compressedParent;
	int code = (int)tree[idx];

	// Compute scalar value for this node
	byte scalar;
	if (code == 0)
		scalar = compressedParent;
	else if (code == 1)
		scalar = (byte)std::min(255.0, parentEstimate + (double)distanceMap[depth]);
	else if (code == 2)
		scalar = (byte)std::max(0.0, parentEstimate - (double)distanceMap[depth]);

	// Recurse on children if not at maximum tree depth
	if (depth < treeDepth) {

			if (depth == leafDepth) {
				compressTreeRecursive(2 * idx + 1, depth + 1, scalar);
				compressTreeRecursive(2 * idx + 2, depth + 1, scalar);
			}
			else {
				compressFromLeaves(2 * idx + 1, depth + 1, scalar, leafDepth);
				compressFromLeaves(2 * idx + 2, depth + 1, scalar, leafDepth);
			}
	}
}