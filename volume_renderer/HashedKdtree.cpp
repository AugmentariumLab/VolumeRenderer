#include "HashedKdtree.h"


inline int64_t HashedKdtree::getCell(int64_t x, int64_t y, int64_t z) {
	return x + (X*y) + (X*Y*z);
}

inline int64_t HashedKdtree::hash(int64_t code) {
	return code & hashMask;
}

inline int64_t leftChildMortonCode(int64_t code) {
	return (code << 1) | 0;
}

inline int64_t rightChildMortonCode(int64_t code) {
	return (code << 1) | 1;
}

void HashedKdtree::build() {

	// Compute height of tree from dimensions of input volume data
	int numXsplits = (int)(log(X) / log(2));
	int numYsplits = (int)(log(Y) / log(2));
	int numZsplits = (int)(log(Z) / log(2));
	treeDepth = numXsplits + numYsplits + numZsplits; // total number of splits (not counting root)
	origTreeDepth = treeDepth;

	// Allocate space for distance map construction
	distanceMap.resize(treeDepth + 1, 0); // treeDepth + root
	distanceSums.resize(treeDepth + 1, 0.0); // treeDepth + root
	distanceCounts.resize(treeDepth + 1, 0.0); // treeDepth + root

	// Allocate space for hashed tree
	numNodes = (int64_t)pow(2, treeDepth + 0); // This has to be high enough to avoid hashing collisions
	hashMask = numNodes - 1;

	temp.resize(numNodes);
	visited.resize(numNodes);

	treeData.resize(numNodes);
	treeStructure.resize(numNodes);
	treeDataCollisions.resize(numNodes);
	treeStructureCollisions.resize(numNodes);
	tempCollisions.resize(numNodes);
	//treeDataCollisions.resize(numNodes / 2);
	//treeStructureCollisions.resize(numNodes / 2);
	//tempCollisions.resize(numNodes / 2);

	// Print tree info
	std::cout << "Tree depth: " << treeDepth << std::endl;
	std::cout << "Number of nodes: " << numNodes << std::endl;
	std::cout << "Tree Size Minimum: " << 2 * ((numNodes * (4.0 / 8.0)) / 1e9) << " GB" << std::endl;

	// Call recursive build function starting from root
	buildRecursive(1, 0, rootMin, rootMax, 0);

	// Compute Distance Map
	for (int depth = 0; depth < treeDepth + 1; depth++) {
		distanceMap[depth] = (byte)(distanceSums[depth] / distanceCounts[depth]);
	}
	std::vector<double>().swap(distanceSums);
	std::vector<double>().swap(distanceCounts);

	std::cout << "BEFORE COMPRESSION COLLISIONS: " << collisionNodes.size() << std::endl;

	// 2nd pass of compression
	compressTreeRecursive(1, 0, 0);

	std::cout << "AFTER COMPRESSION COLLISIONS: " << collisionNodes.size() << std::endl;
	numCollisions = collisionNodes.size();

	// Print Distance Map
	for (int depth = 0; depth < treeDepth + 1; depth++) {
		std::cout << (int)distanceMap[depth] << std::endl;
	}

	// Print Non-empty nodes
	int64_t nonEmpties = 0;
	for (int64_t node = 0; node < numNodes; node++) {
		if (treeData[node] != 3)
			nonEmpties++;
	}
	std::cout << "NON-EMPTY NODES: " << nonEmpties << std::endl;
	std::cout << "EMPTY NODES: " << numNodes - nonEmpties << std::endl;

	// delete temporary trees
	std::vector<byte>().swap(temp);
	std::vector<byte>().swap(tempCollisions);
	std::vector<int64_t>().swap(visited);

	treeDataCollisions.resize(numCollisions);
	treeStructureCollisions.resize(numCollisions);
}

void HashedKdtree::buildRecursive(int64_t mcode, int depth, Point3i minBound,
	Point3i maxBound, byte parentEstimate) {

	int64_t key = hash(mcode);
	bool isCollision = (treeData[key] == 3);
	byte scalarEstimate;

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

	// Create a leaf node if all scalars are uniform within bounding box
	bool isLeaf = false;
	if (minScalar == maxScalar) {
		isLeaf = true;
	}

	// Compute midrange scalar value of all volume cells
	byte midrange = (byte)((maxScalar + minScalar) / 2.0);

	// Determine if node is a collision with other nodes according to morton key
	if (!isCollision) {
		if (visited[key] == 0) {
			// this key has not be visited, NOT a collision
			visited[key] = mcode;
		}
		else if (visited[key] != mcode) {
			// this key has been visited by another tree node (COLLISION FOUND!)
			
			if (collisionNodes.count(visited[key]) == 0) {
				int64_t prevNodeCollisionIdx = lastCollisionIdx++;
				collisionNodes.insert({ { visited[key], prevNodeCollisionIdx } });
				tempCollisions[prevNodeCollisionIdx] = temp[key];
				treeStructureCollisions[prevNodeCollisionIdx] = (int)treeStructure[key];
				treeDataCollisions[prevNodeCollisionIdx] = (int)treeData[key];
			}

			if (collisionNodes.count(mcode) == 0) {
				int64_t currentNodeCollisionIdx = lastCollisionIdx++;
				collisionNodes.insert({ { mcode, currentNodeCollisionIdx } });
			}
			treeData[key] = 3;
			isCollision = true;
		}
	}
	

	if (isCollision) {
		if (collisionNodes.count(mcode) == 0) {
			int64_t currentNodeCollisionIdx = lastCollisionIdx++;
			collisionNodes.insert({ { mcode, currentNodeCollisionIdx } });
		}
		key = collisionNodes[mcode]; // collisionIdx
		tempCollisions[key] = midrange;
		scalarEstimate = encodeNode(mcode, depth, parentEstimate, midrange, false, key);
	}
	else {
		temp[key] = midrange;
		scalarEstimate = encodeNode(mcode, depth, parentEstimate, midrange, false);
	}
	

	// Recurse on children if not at maximum tree depth
	if (depth < treeDepth) {

		if (!isLeaf)
			if (isCollision)
				treeStructureCollisions[collisionNodes[mcode]] = 3;
			else
				treeStructure[key] = 3;

		Point3i extent = maxBound - minBound;
		int64_t numCells = extent.prod();
		
		// Append bit to current node's morton code to get morton code of children
		int64_t leftCode = leftChildMortonCode(mcode);
		int64_t rightCode = rightChildMortonCode(mcode);

		// Do NOT split bounding box if only 1 cell is enclosed
		if (numCells == 1) {
			buildRecursive(leftCode, depth + 1, minBound, maxBound, scalarEstimate); // left
			buildRecursive(rightCode, depth + 1, minBound, maxBound, scalarEstimate); // right
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
			buildRecursive(leftCode, depth + 1, minBound, maxBound, scalarEstimate);

			// right child
			minBound[splitDim] = thisMid;
			maxBound[splitDim] = thisMax;
			buildRecursive(rightCode, depth + 1, minBound, maxBound, scalarEstimate);

		}
	}
}

void HashedKdtree::levelCut(int cutDepth, std::vector<byte> &outData) {

	output = &outData;
	queryDepth = cutDepth;

	// Resize output vector
	output->resize(X*Y*Z);

	// Call recursive level cut function from root
	levelCutRecursive(1, 0, rootMin, rootMax, 0);
}

int HashedKdtree::measureMaxError() {
	int maxError = 0;
	for (int64_t i = 0; i < X*Y*Z; i++) {
		maxError = std::max((int)abs((double)(*output)[i] - (double)(*data)[i]), maxError);
	}
	return maxError;
}

double HashedKdtree::measureMeanError() {
	double sumError = 0.0;
	int64_t count = X*Y*Z;
	for (int64_t i = 0; i < count; i++) {
		sumError += abs((double)(*output)[i] - (double)(*data)[i]);
	}
	return sumError / (double)count;
}

void HashedKdtree::queryError(std::vector<byte> &outData) {
	output2 = &outData;
	int64_t count = X*Y*Z;
	output2->resize(count);
	for (int64_t i = 0; i < count; i++) {
		(*output2)[i] = (byte)abs((double)(*output)[i] - (double)(*data)[i]);
	}
}

void HashedKdtree::levelCutRecursive(int64_t mcode, int depth, Point3i minBound, Point3i maxBound, byte compressedParent) {

	double parentEstimate = (double)compressedParent;
	int64_t key = hash(mcode);
	int code = (int)treeData[key];

	// Compute scalar value for this node
	byte scalar;
	bool isCollision = (code == 3);
	if (isCollision) {
		key = collisionNodes[mcode];
		code = treeDataCollisions[key];
	}
		
	if (code == 0)
		scalar = compressedParent;
	else if (code == 1)
		scalar = (byte)std::min(255.0, parentEstimate + (double)distanceMap[depth]);
	else if (code == 2)
		scalar = (byte)std::max(0.0, parentEstimate - (double)distanceMap[depth]);


	// Populate output vectors if at cut depth or a leaf
	int children;
	if (isCollision)
		children = treeStructureCollisions[key];
	else
		children = treeStructure[key];

	if (depth == queryDepth || children == 0) {
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
		Point3i extent = maxBound - minBound;
		int64_t numCells = extent.prod();

		// Append bit to current node's morton code to get morton code of children
		int64_t leftCode = leftChildMortonCode(mcode);
		int64_t rightCode = rightChildMortonCode(mcode);

		// Do NOT split bounding box if only 1 cell is enclosed
		if (numCells == 1) {
			if (children == 3 || children == 1)
				levelCutRecursive(leftCode, depth + 1, minBound, maxBound, scalar); // left
			if (children == 3 || children == 2)
				levelCutRecursive(rightCode, depth + 1, minBound, maxBound, scalar); // right
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
			if (children == 3 || children == 1) {
				maxBound[splitDim] = thisMid;
				levelCutRecursive(leftCode, depth + 1, minBound, maxBound, scalar);
			}
			// right child
			if (children == 3 || children == 2) {
				minBound[splitDim] = thisMid;
				maxBound[splitDim] = thisMax;
				levelCutRecursive(rightCode, depth + 1, minBound, maxBound, scalar);
			}
		}
}

byte HashedKdtree::encodeNode(int64_t mcode, int depth, byte compressedParent, byte trueNode, bool useMap, int collisionIdx) {

	// For easy reference
	double parentEstimate = (double)compressedParent;
	double nodeTruth = (double)trueNode;
	//double nodeTruth = (double)temp[hash(mcode)];
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
		if (useMap)
			if (collisionIdx > -1)
				treeDataCollisions[collisionIdx] = 0;
			else
				treeData[hash(mcode)] = 0;
		return (byte)noneEstimate;
	}

	// If adding produces minimum error, set node to 1 & return corresponding estimate
	if (minError == addError) {
		if (!useMap) {
			distanceSums[depth] += addError;
			distanceCounts[depth]++;
		}
		else
			if (collisionIdx > -1)
				treeDataCollisions[collisionIdx] = 1;
			else
				treeData[hash(mcode)] = 1;
		return (byte)addEstimate;
	}

	// If subtracting produces minimum error, set node to 2 & return corresponding estimate
	if (minError == subError) {
		if (!useMap) {
			distanceSums[depth] += subError;
			distanceCounts[depth]++;
		} else
			if (collisionIdx > -1)
				treeDataCollisions[collisionIdx] = 2;
			else
				treeData[hash(mcode)] = 2;
		return (byte)subEstimate;
	}
}

void HashedKdtree::compressTreeRecursive(int64_t mcode, int depth, byte compressedParent, int trueNodeOverride) {
	
	int64_t key = hash(mcode);
	byte compressedScalar, trueScalar;
	int children = 0;
	bool isCollision = (treeData[key] == 3);

	// Determine if node is a collision with other nodes according to morton key
	if (!isCollision) {
		if (visited[key] == 0) {
			// this key has not be visited, NOT a collision
			visited[key] = mcode;
		}
		else if (visited[key] != mcode) {
			// this key has been visited by another tree node (COLLISION FOUND!)
			if (collisionNodes.count(visited[key]) == 0) {
				int64_t prevNodeCollisionIdx = lastCollisionIdx++;
				collisionNodes.insert({ { visited[key], prevNodeCollisionIdx } });
				tempCollisions[prevNodeCollisionIdx] = temp[key];
				treeStructureCollisions[prevNodeCollisionIdx] = (int)treeStructure[key];
				treeDataCollisions[prevNodeCollisionIdx] = (int)treeData[key];
			}
			if (collisionNodes.count(mcode) == 0) {
				int64_t currentNodeCollisionIdx = lastCollisionIdx++;
				collisionNodes.insert({ { mcode, currentNodeCollisionIdx } });
			}
			treeData[key] = 3;
			isCollision = true;
		}
	}


	if (isCollision) {
		if (collisionNodes.count(mcode) == 0) {
			int64_t currentNodeCollisionIdx = lastCollisionIdx++;
			collisionNodes.insert({ { mcode, currentNodeCollisionIdx } });
		}
		key = collisionNodes[mcode]; // collisionIdx
		children = treeStructureCollisions[key];
		trueScalar = trueNodeOverride != -1 ? (byte)trueNodeOverride : tempCollisions[key];
		compressedScalar = encodeNode(mcode, depth, compressedParent, trueScalar, true, key);
	}
	else {
		children = treeStructure[key];
		trueScalar = trueNodeOverride != -1 ? (byte)trueNodeOverride : temp[key];
		compressedScalar = encodeNode(mcode, depth, compressedParent, trueScalar, true);
	}

	if (children == 0) {

		int leafError = abs((int)compressedScalar - (int)trueScalar);

		// Split into 2 children if above original depth & any error present
		if (leafError > 0 && depth < origTreeDepth) {
			if (isCollision)
				children = treeStructureCollisions[key] = 3;
			else
				children = treeStructure[key] = 3;
		}

		// If leaf error is too high & there is room to grow the branch, grow the branch
		else if (leafError > tolerance && (depth < treeDepth || maxAddLevels > 0)) {

			//std::cout << (int)scalar << " " << (int)trueNode << std::endl; // DEBUG

			// randomly choose a child to grow tree branch
			// TODO: Intentionally choose child to avoid collisions!
			std::shuffle(childVec.begin(), childVec.end(), gen);
			if (isCollision)
				children = treeStructureCollisions[key] = childVec[0];
			else
				children = treeStructure[key] = childVec[0];
			//treeStructure[key] = childVec[0];
			//children = treeStructure[key];
			//std::cout << children << std::endl; // DEBUG

			// carry over true node value to child recursions
			trueNodeOverride = (int)trueScalar;
			//std::cout << trueNodeOverride << std::endl; // DEBUG

			// Add level to tree if reached max depth
			if (depth == treeDepth && maxAddLevels > 0) {
				// add level to tree & distance map
				treeDepth++;
				distanceMap.resize(treeDepth + 1, addLevelDistance);
				addLevelDistance /= 2;
				maxAddLevels--;
				//std::cout << (int)distanceMap[treeDepth] << std::endl; // DEBUG
			}
		}
		else {
			return;
		}
	}
		
	// Recurse on left child
	if (children == 3 || children == 1)
		compressTreeRecursive(leftChildMortonCode(mcode), depth + 1, compressedScalar, trueNodeOverride);
	// Recurse on right child
	if (children == 3 || children == 2)
		compressTreeRecursive(rightChildMortonCode(mcode), depth + 1, compressedScalar, trueNodeOverride);
}

void HashedKdtree::save(std::string filename) {
	// Get size of  KDOctree 
	int64_t treeSize = treeData.bytes();

	// write keys and values of unordered map
	std::vector<int64_t> keys;
	keys.reserve(numCollisions);
	std::vector<int32_t> vals;
	vals.reserve(numCollisions);
	for (auto kv : collisionNodes) {
		keys.push_back(kv.first);
		vals.push_back(kv.second);
	}
	int64_t keysSize = numCollisions * sizeof(int64_t);
	int64_t valsSize = numCollisions * sizeof(int32_t);

	// Write out the data
	std::ofstream out(filename, std::ios::out | std::ios::binary);

	out.write(reinterpret_cast<char *>(&rootMin), sizeof(Point3i));
	out.write(reinterpret_cast<char *>(&rootMax), sizeof(Point3i));

	out.write(reinterpret_cast<char *>(&treeDepth), sizeof(int));
	out.write(reinterpret_cast<char *>(&X), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Z), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&hashMask), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&numCollisions), sizeof(int64_t)); // number of collisions

	out.write(reinterpret_cast<char *>(&distanceMap[0]), treeDepth + 1);
	out.write(reinterpret_cast<char *>(&treeData.bits[0]), treeSize);
	out.write(reinterpret_cast<char *>(&treeStructure.bits[0]), treeSize);
	out.write(reinterpret_cast<char *>(&treeDataCollisions.bits[0]), (numCollisions + 3) / 4);
	out.write(reinterpret_cast<char *>(&treeStructureCollisions.bits[0]), (numCollisions + 3) / 4);

	out.write(reinterpret_cast<char *>(&keys[0]), keysSize);
	out.write(reinterpret_cast<char *>(&vals[0]), valsSize);

	out.close();

	// Print info
	int64_t dataSize = keysSize + valsSize + (2 * treeSize) + (2 * numCollisions) + 
		(2 * sizeof(Point3i)) + sizeof(int) + treeDepth + 1 + (5 * sizeof(int64_t));
	std::cout << "\nNew file saved: " << filename << std::endl;
	std::cout << "File size: " << (double)dataSize / 1e9 << " GB" << std::endl;
}

void HashedKdtree::open(std::string filename) {

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
	is.read(reinterpret_cast<char *>(&hashMask), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&numCollisions), sizeof(int64_t)); // number of collisions

	int64_t collisionsSize = (numCollisions + 3) / 4;
	int64_t treeSize = fileSize - ( 
		(2 * sizeof(Point3i)) + 
		sizeof(int) + 
		(treeDepth + 1) + 
		(5 * sizeof(int64_t)) + 
		(collisionsSize * 2) + 
		(numCollisions * sizeof(int64_t)) + 
		(numCollisions * sizeof(int32_t)) 
		);
	treeSize /= 2;
	std::cout << treeSize << std::endl;

	// Allocate memory for distanceMap & tree
	distanceMap.resize(treeDepth + 1);
	treeData.bits.resize(treeSize);
	treeStructure.bits.resize(treeSize);
	treeDataCollisions.resize(numCollisions);
	treeStructureCollisions.resize(numCollisions);
	std::vector<int64_t> keys (numCollisions);
	std::vector<int32_t> vals (numCollisions);

	is.read(reinterpret_cast<char *>(&distanceMap[0]), treeDepth + 1);
	is.read(reinterpret_cast<char *>(&treeData.bits[0]), treeSize);
	is.read(reinterpret_cast<char *>(&treeStructure.bits[0]), treeSize);
	is.read(reinterpret_cast<char *>(&treeDataCollisions.bits[0]), (numCollisions + 3) / 4);
	is.read(reinterpret_cast<char *>(&treeStructureCollisions.bits[0]), (numCollisions + 3) / 4);

	is.read(reinterpret_cast<char *>(&keys[0]), numCollisions * sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&vals[0]), numCollisions * sizeof(int32_t));

	// populate map
	for (int i = 0; i < numCollisions; i++) {
		collisionNodes.emplace(keys[i], vals[i]);
	}

	// Close the file
	is.close();

}
