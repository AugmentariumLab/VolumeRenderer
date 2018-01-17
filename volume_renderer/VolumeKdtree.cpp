#include "VolumeKdtree.h"


inline int64_t VolumeKdtree::getCell(int64_t x, int64_t y, int64_t z) {
	return x + (X*y) + (X*Y*z);
}


void VolumeKdtree::setErrorTolerance(int errorTolerance) {
	tolerance = errorTolerance;
}

void VolumeKdtree::setMaxEpochs(int epochs) {
	maxEpochs = epochs;
}

void VolumeKdtree::build(bool useThreads) {

	parallelism = useThreads;

	// Define additional distance map levels
	int maxAddLevels = 7;
	int addLevelDistance[] = { 64, 32, 16, 8, 4, 2, 1 };

	// Analyze input data dimensions to define size of Kd-tree ////////////////////////
	int numXsplits = (int)(log(X) / log(2));
	int numYsplits = (int)(log(Y) / log(2));
	int numZsplits = (int)(log(Z) / log(2));
	origTreeDepth = numXsplits + numYsplits + numZsplits; // total number of splits (not counting root)
	maxTreeDepth = origTreeDepth + maxAddLevels; 
	distanceMap.resize(maxTreeDepth + 1, 0); // treeDepth + root
	distanceSums.resize(maxTreeDepth + 1, 0.0); // treeDepth + root
	distanceCounts.resize(maxTreeDepth + 1, 0.0); // treeDepth + root
	numOrigNodes = (int64_t)pow(2, origTreeDepth + 1) - 1;
	numMaxNodes = numOrigNodes + (int64_t)(pow(2, origTreeDepth)* (maxAddLevels));
	temp.resize(numOrigNodes);

	// Print data/tree info
	std::cout << "Original Dataset size: " << double(data->size() * sizeof((*data)[0])) / 1e9 << " GB" << std::endl;
	std::cout << "Maximum Tree depth: " << maxTreeDepth << std::endl;
	std::cout << "Maximum Tree size: " << (numMaxNodes * (2.0 / 8.0)) / 1e9 << " GB" << std::endl;


	// BUILD //////////////////////////////////////////////////////////////////////////
	// Call recursive build function starting from root
	// Populates 'temp' tree array (byte array)
	DebugTimer::Begin(1, "BUILD");
	buildRecursive(0, 0, rootMin, rootMax);

	// Deallocate original dataset
	data->clear();
	data->shrink_to_fit();
	DebugTimer::End("BUILD");

	// COMPRESS /////////////////////////////////////////////////////////////////////
	// Define distance map & compress 'temp' byte array into 2-bit 
	// 'tree' array using progressive scheme
	DebugTimer::Begin(1, "COMPRESS");
	tree.resize(numOrigNodes); // 2-bit array elements
	firstOrigLeaf = (int64_t)pow(2, origTreeDepth) - 1; // Index of 1st leaf in breadth-first arrays
	compressGradientDescent();

	// Reduce temporary tree to contain only leaves
	temp.erase(temp.begin(), temp.begin() + firstOrigLeaf);
	temp.shrink_to_fit();
	int64_t numLeaves = temp.size();
	DebugTimer::End("COMPRESS");
	
	

	// Max Error
	int maxError = 0;
	for (int64_t i = 0; i < numLeaves; i++) {
		maxError = std::max((int)abs((double)temp[i] - (double)recon[i]), maxError);
	}
	std::cout << "\nBefore branch growth, Maximum Error: " << maxError << std::endl;

	double sumError = 0.0;
	double sumErrorL2 = 0.0;
	for (int64_t i = 0; i < numLeaves; i++) {
		sumError += abs((double)temp[i] - (double)recon[i]);
		sumErrorL2 += pow( ((double)temp[i] - (double)recon[i]), 2.0);
	}
	std::cout << "Before branch growth, Mean Error L1: " << sumError / (double)numLeaves << "  Mean Error L2: " << sumErrorL2 / (double)numLeaves << std::endl;

	// PRUNE /////////////////////////////////////////////////////////////////////////
	// Prune Tree by converting leaf nodes to '3' code
	// 3rd PASS - uses stack
	DebugTimer::Begin(1, "PRUNE");
	pruneTreeRecursive(0);
	DebugTimer::End("PRUNE");

	// Fix additional levels
	int addIdx = 0;
	for (int depth = origTreeDepth + 1; depth < maxTreeDepth + 1; depth++) {
		distanceMap[depth] = (byte)addLevelDistance[addIdx++];
	}

	// Print distance map - DEBUG ///
	/*
	std::cout << "DISTANCE MAP: " << std::endl;
	for (int depth = 0; depth < maxTreeDepth + 1; depth++) {
		std::cout << (int)distanceMap[depth] << std::endl;
	}
	*/

	// CONVERT //////////////////////////////////////////////////////////////////////
	// Convert breadth-first 'tree' array to unbalanced depth-first 'tree' array
	// for better compression. Grows branches if leaf over error tolerance.
	// 4th PASS - uses stack
	DebugTimer::Begin(1, "CONVERT");
	convertToPreorder();
	DebugTimer::Begin(1, "CONVERT");

	// Max Error
	maxError = 0;
	for (int64_t i = 0; i < numLeaves; i++) {
		maxError = std::max((int)abs((double)temp[i] - (double)recon[i]), maxError);
	}
	std::cout << "\nMaximum Error: " << maxError << std::endl;

	// Mean Error
	sumErrorL2 = 0.0;
	sumError = 0.0;
	for (int64_t i = 0; i < numLeaves; i++) {
		sumError += abs((double)temp[i] - (double)recon[i]);
		sumErrorL2 += pow(((double)temp[i] - (double)recon[i]), 2.0);
	}
	std::cout << "Mean L1 Error: " << sumError / (double)numLeaves << "  Mean Error L2: " << sumErrorL2 / (double)numLeaves << std::endl;
	
	// Deallocate temporary trees
	//temp.clear();
	//temp.shrink_to_fit();
	//recon.clear();
	//recon.shrink_to_fit();

	// Print compressed tree info
	std::cout << "\nNumber of active nodes: " << numActiveNodes << std::endl;
	std::cout << "Compressed Tree Size : " << (numActiveNodes * (2.0 / 8.0)) / 1e9 << " GB" << std::endl;
}


MinMax VolumeKdtree::buildRecursive(int64_t idx, int depth, Point3i minBound, Point3i maxBound) {

	MinMax leftScalarRange, rightScalarRange;
	double thisMinScalar, thisMaxScalar;

	// Recurse on children if not at maximum tree depth
	if (depth < origTreeDepth) {

		int splitDim = depth % MAX_DIM;
		Point3i extent = maxBound - minBound;
		int64_t numCells = extent.prod();

		// Choose split dimension to acheive full resolution
		int i = 0;
		while (numCells > 1 && extent[splitDim] == 1) {
			splitDim = (depth + ++i) % MAX_DIM;
		}

		// Split Bounding box 
		int64_t thisMid = (minBound[splitDim] + maxBound[splitDim]) / 2;
		int64_t thisMax = maxBound[splitDim];

		if (parallelism) {
			maxBound[splitDim] = thisMid;
			Point3i leftMin = minBound;
			Point3i leftMax = maxBound;

			minBound[splitDim] = thisMid;
			maxBound[splitDim] = thisMax;
			Point3i rightMin = minBound;
			Point3i rightMax = maxBound;

			parallel_invoke(
				[this, idx, depth, minBound, leftMin, leftMax, &leftScalarRange] { leftScalarRange = this->buildRecursive(2 * idx + 1, depth + 1, leftMin, leftMax); },
				[this, idx, depth, minBound, rightMin, rightMax, &rightScalarRange] { rightScalarRange = this->buildRecursive(2 * idx + 2, depth + 1, rightMin, rightMax);  }
			);
		}
		else {
			// left child
			maxBound[splitDim] = thisMid;
			leftScalarRange = buildRecursive(2 * idx + 1, depth + 1, minBound, maxBound);

			// right child
			minBound[splitDim] = thisMid;
			maxBound[splitDim] = thisMax;
			rightScalarRange = buildRecursive(2 * idx + 2, depth + 1, minBound, maxBound);
		}

		thisMinScalar = std::min((double)leftScalarRange.min, (double)rightScalarRange.min);
		thisMaxScalar = std::max((double)leftScalarRange.max, (double)rightScalarRange.max);
	} 
	else if (depth == origTreeDepth)
		thisMaxScalar = thisMinScalar = (*data)[getCell(minBound[0], minBound[1], minBound[2])];

	// Populate true tree array
	temp[idx] = (byte)((thisMaxScalar + thisMinScalar) / 2.0); //midrange

	return MinMax((byte)thisMinScalar, (byte)thisMaxScalar);
}

/**
Find optimal distance map values for each tree level using gradient descent.
*/
void VolumeKdtree::compressGradientDescent() {

	// Constant variables -- gradient descent knobs
	const double gamma = 1.25; // step size multiplier
	const double h = 1.0; // central difference interval
	const double maxAbsStepSize = 4.0;

	// Declare processing variables

	// These variables are updated @ each tree depth
	//int64_t nodeIdx, 
	int64_t numNodes;
	int64_t startingNodeIdx = 0, endingNodeIdx = 0, parentStartingNodeIdx = 0;
	std::vector<byte> reconstructedParents; // reconstructed node @ depth - 1
	std::vector<byte> reconPreviousEpoch; // save recon from previous epoch in case of revert

	// These variables are updated @ each gradient descent epoch
	int epoch;
	double distanceSum, distanceCount;
	double currentDistance, currentError, currentDF, currentStepSize;
	double previousStepSize, previousDistance, previousDF, previousError;
	//byte reconstructedParent, reconstructedNode, trueNode;

	for (int depth = 0; depth < origTreeDepth + 1; depth++) {
		
		// Initialize processing variables
		epoch = 0;
		distanceSum = distanceCount = 0.0;
		numNodes = (int64_t)pow(2, depth);
		endingNodeIdx = startingNodeIdx + numNodes;
		recon.resize(numNodes);

		// Find starting distance for gradient descent ////////////////////////////////////

		// Loop through nodes in level order to populate distanceSum & distanceCount
		//if (parallelism) {
		if (false) { // NOT WORKING!
			combinable<double> sum([]() { return 0; });
			combinable<double> count([]() { return 0; });
			parallel_for(int64_t(startingNodeIdx), endingNodeIdx, [this, reconstructedParents,
				parentStartingNodeIdx, &sum, &count](int64_t nodeIdx) {
				byte reconstructedParent = nodeIdx == 0 ? 0 : reconstructedParents[((nodeIdx - 1) / 2) - parentStartingNodeIdx];
				encodeNodeEstimate(nodeIdx, reconstructedParent, &sum.local(), &count.local());
			});
			distanceSum = sum.combine(plus<double>());
			distanceCount = count.combine(plus<double>());
		} 
		else {
			for (int64_t nodeIdx = startingNodeIdx; nodeIdx < endingNodeIdx; nodeIdx++) {
				byte reconstructedParent = nodeIdx == 0 ? 0 : reconstructedParents[((nodeIdx - 1) / 2) - parentStartingNodeIdx];
				encodeNodeEstimate(nodeIdx, reconstructedParent, &distanceSum, &distanceCount);
			}
		}
		


		// Calculate current distance value from distanceSum & distanceCount
		if (distanceCount > 0)
			currentDistance = round(distanceSum / distanceCount);
		else
			currentDistance = 0.0;
		
		// End find starting distance block //////////////////////////////////////////////


		// Perform gradient descent from starting distance //////////////////////////////
		previousDistance = 0.0;
		previousStepSize = 255.0;
		previousError = 65025.0;
		while (epoch < maxEpochs && abs(previousStepSize) >= 0.5) {

			// Update variables
			if (epoch != 0) {
				previousDistance = currentDistance;
				previousError = currentError;
				previousDF = currentDF;
				previousStepSize = currentStepSize;

				currentDistance = round(std::min(255.0, std::max(0.0, previousDistance + previousStepSize)));

				// do NOT re-evaluate error & df if currentDistance == previousDistance
				if (currentDistance == previousDistance)
					break;
			}
			
			// Evaluate error for current distance
			//if (parallelism) {
			if (false) {
				combinable<double> error([]() { return 0; });
				parallel_for(int64_t(startingNodeIdx), endingNodeIdx, [this, reconstructedParents,
					parentStartingNodeIdx, currentDistance, &error, startingNodeIdx](int64_t nodeIdx) {
					byte reconstructedParent = nodeIdx == 0 ? 0 : reconstructedParents[((nodeIdx - 1) / 2) - parentStartingNodeIdx];
					double L1Error;
					byte reconstructedNode = encodeNode(nodeIdx, reconstructedParent, (byte)currentDistance, true, -1.0, &L1Error);
					//currentError += pow(L1Error, 2.0);
					error.local() += pow(L1Error, 2.0);
					recon[nodeIdx - startingNodeIdx] = reconstructedNode; 
				});
				currentError = error.combine(plus<double>());
			}
			else {
				for (int64_t nodeIdx = startingNodeIdx; nodeIdx < endingNodeIdx; nodeIdx++) {
					byte reconstructedParent = nodeIdx == 0 ? 0 : reconstructedParents[((nodeIdx - 1) / 2) - parentStartingNodeIdx];
					double L1Error;
					byte reconstructedNode = encodeNode(nodeIdx, reconstructedParent, (byte)currentDistance, true, -1.0, &L1Error);
					currentError += pow(L1Error, 2.0);
					recon[nodeIdx - startingNodeIdx] = reconstructedNode;
				}
			}
			currentError /= numNodes;
			

			// Break gradient descent loop if found zero-error distance map value
			if (currentError < 1.0)
				break;

			// Revert back if currentError  > previousError
			if (epoch != 0 && currentError > previousError) {
				currentError = previousError;
				currentDistance = previousDistance;
				currentDF = previousDF;
				currentStepSize = previousStepSize / 2.0;
				recon.swap(reconPreviousEpoch);
				epoch++;
				continue;
			}

			// Evaluate DF for current distance using central difference fomula
			byte estimateDistances[2] = { (byte)std::max(0.0, currentDistance - h), (byte)std::min(255.0, currentDistance + h) };
			double estimateErrors[2] = { 0, 0 };
			if (parallelism) {
			//if (false) {
				parallel_for(int(0), 2, [this, &estimateDistances, &estimateErrors, startingNodeIdx,
					endingNodeIdx, reconstructedParents, parentStartingNodeIdx, numNodes](int i) {
					for (int64_t nodeIdx = startingNodeIdx; nodeIdx < endingNodeIdx; nodeIdx++) {
						byte reconstructedParent = nodeIdx == 0 ? 0 : reconstructedParents[((nodeIdx - 1) / 2) - parentStartingNodeIdx];
						double L1error;
						encodeNode(nodeIdx, reconstructedParent, estimateDistances[i], false, -1.0, &L1error);
						estimateErrors[i] += pow(L1error, 2.0);
					}
					estimateErrors[i] /= numNodes;
				});
			}
			else {
				for (int i = 0; i < 2; i++) {
					for (int64_t nodeIdx = startingNodeIdx; nodeIdx < endingNodeIdx; nodeIdx++) {
						byte reconstructedParent = nodeIdx == 0 ? 0 : reconstructedParents[((nodeIdx - 1) / 2) - parentStartingNodeIdx];
						double L1error;
						encodeNode(nodeIdx, reconstructedParent, estimateDistances[i], false, -1.0, &L1error);
						estimateErrors[i] += pow(L1error, 2.0);
					}
					estimateErrors[i] /= numNodes;
				}
			}

			currentDF = (estimateErrors[1] - estimateErrors[0]) / (2.0 * h);
			currentStepSize = std::max(-maxAbsStepSize, std::min(maxAbsStepSize, -gamma * currentDF));
			//std::cout << currentDF << " " << currentError << " " << currentDistance << " " << depth << " " << epoch << std::endl; //DEBUG 
			reconPreviousEpoch = recon;
			epoch++;
		}

		// set distance map to optimal distance
		this->distanceMap[depth] = (byte)currentDistance;
		//std::cout << currentError << std::endl;
		//std::cout << " --> " << (int)distanceMap[depth] << std::endl; // DEBUG

		// Reset variables
		if (depth < origTreeDepth) {
			reconstructedParents.swap(recon);
			recon.clear();
			parentStartingNodeIdx = startingNodeIdx;
			startingNodeIdx = endingNodeIdx;
			reconPreviousEpoch.clear();
			reconPreviousEpoch.shrink_to_fit();
		}
	}

}

int VolumeKdtree::measureMaxError() {
	int maxError = 0;
	for (int64_t i = 0; i < X*Y*Z; i++) {
		maxError = std::max((int)abs((double)(*output)[i] - (double)(*data)[i]), maxError);
	}
	return maxError;
}

double VolumeKdtree::measureMeanError() {
	double sumError = 0.0;
	int64_t count = X*Y*Z;
	for (int64_t i = 0; i < count; i++) {
		sumError += abs((double)(*output)[i] - (double)(*data)[i]);
	}
	return sumError / (double)count;
}


void VolumeKdtree::queryError(std::vector<byte> &outData) {
	output2 = &outData;
	int64_t count = X*Y*Z;
	output2->resize(count);
	for (int64_t i = 0; i < count; i++) {
		(*output2)[i] = (byte)abs((double)(*output)[i] - (double)(*data)[i]);
	}
}



byte VolumeKdtree::encodeNodeEstimate(int64_t idx, byte compressedParent, double * estimateSum, double * estimateCount) {
	
	double nodeTruth = (double)temp[idx];
	double parentEstimate = (double)compressedParent;
	double parentDistance = abs(parentEstimate - nodeTruth);
	double masterDistance = (*estimateSum + parentDistance) / (*estimateCount + 1.0);

	// Calculate error of *doing nothing* 
	double noneEstimate = parentEstimate;
	double noneError = parentDistance;

	// Calculate error of *adding* the master distance
	double addEstimate = std::min(255.0, parentEstimate + masterDistance);
	double addError = abs(addEstimate - nodeTruth);

	// Calculate error of *subtracting* the master difference
	double subEstimate = std::max(0.0, parentEstimate - masterDistance);
	double subError = abs(subEstimate - nodeTruth);

	// Determine minimum error scenario
	double minError = std::min(subError, std::min(noneError, addError));

	// If doing nothing produces minimum error, set node to 0 & return corresponding estimate
	if (minError == noneError) {
		return (byte)noneEstimate;
	}

	// If adding produces minimum error, return corresponding estimate
	if (minError == addError) {
		*estimateSum += parentDistance;
		*estimateCount += 1.0;
		return (byte)addEstimate;
	}

	// If subtracting produces minimum error, return corresponding estimate
	if (minError == subError) {
		*estimateSum += parentDistance;
		*estimateCount += 1.0;
		return (byte)subEstimate;
	}
}

byte VolumeKdtree::encodeNode(int64_t idx, byte compressedParent, byte distanceVal, bool fillTree, double nodeTruth, double * returnError) {

	if (nodeTruth == -1.0)
		nodeTruth = (double)temp[idx];
	double parentEstimate = (double)compressedParent;
	double parentDistance = abs(parentEstimate - nodeTruth);
	double masterDistance = (double)distanceVal;

	// Calculate error of *doing nothing* 
	double noneEstimate = parentEstimate;
	double noneError = parentDistance;

	// Calculate error of *adding* the master distance
	double addEstimate = std::min(255.0, parentEstimate + masterDistance);
	double addError = abs(addEstimate - nodeTruth);

	// Calculate error of *subtracting* the master difference
	double subEstimate = std::max(0.0, parentEstimate - masterDistance);
	double subError = abs(subEstimate - nodeTruth);

	// Determine minimum error scenario
	double minError = std::min(subError, std::min(noneError, addError));
	if (returnError)
		*returnError = minError;

	// If doing nothing produces minimum error, set node to 0 & return corresponding estimate
	if (minError == noneError) {
		if (fillTree)
			tree[idx] = 0;
		return (byte)noneEstimate;
	}

	// If adding produces minimum error, return corresponding estimate
	if (minError == addError) {
		if (fillTree)
			tree[idx] = 1;
		return (byte)addEstimate;
	}

	// If subtracting produces minimum error, return corresponding estimate
	if (minError == subError) {
		if (fillTree)
			tree[idx] = 2;
		return (byte)subEstimate;
	}
}
/*
void VolumeKdtree::compressTreeRecursive(int64_t idx, int depth, byte compressedParent) {

	// Find scalar value & code for current node
	byte scalar = encodeNode(idx, depth, compressedParent, true);
	if (depth == origTreeDepth)
		recon[idx - firstOrigLeaf] = scalar;

	// Recurse on children
	//if (depth < treeDepth) {
	else if (depth < origTreeDepth) {
		parallel_invoke(
			[this, idx, depth, scalar] { this->compressTreeRecursive(2 * idx + 1, depth + 1, scalar); },
			[this, idx, depth, scalar] { this->compressTreeRecursive(2 * idx + 2, depth + 1, scalar); }
		);
		//compressTreeRecursive(2*idx + 1, depth + 1, scalar);
		//compressTreeRecursive(2*idx + 2, depth + 1, scalar);
	}
}
*/
void VolumeKdtree::save(std::string filename) {
	// Get size of  KD-tree 
	int64_t treeSize = tree.bytes();

	// Exit if tree id empty
	if (treeSize == 0) {
		std::cout << "ERROR! No tree to save." << std::endl;
		std::cin.ignore();
		return;
	}

	// Write out the data
	std::ofstream out(filename, std::ios::out | std::ios::binary);

	out.write(reinterpret_cast<char *>(&rootMin), sizeof(Point3i));
	out.write(reinterpret_cast<char *>(&rootMax), sizeof(Point3i));
	out.write(reinterpret_cast<char *>(&maxTreeDepth), sizeof(int));
	out.write(reinterpret_cast<char *>(&origTreeDepth), sizeof(int));
	out.write(reinterpret_cast<char *>(&X), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&Z), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&numActiveNodes), sizeof(int64_t));
	out.write(reinterpret_cast<char *>(&distanceMap[0]), maxTreeDepth + 1);
	out.write(reinterpret_cast<char *>(&tree.bits[0]), treeSize);

	out.close();

	// Print info
	int64_t dataSize = treeSize + (2 * sizeof(Point3i)) + (2 * sizeof(int)) + maxTreeDepth + 1 + (4 * sizeof(int64_t));
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

	is.read(reinterpret_cast<char *>(&rootMin), sizeof(Point3i));
	is.read(reinterpret_cast<char *>(&rootMax), sizeof(Point3i));
	is.read(reinterpret_cast<char *>(&maxTreeDepth), sizeof(int));
	is.read(reinterpret_cast<char *>(&origTreeDepth), sizeof(int));
	is.read(reinterpret_cast<char *>(&X), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&Y), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&Z), sizeof(int64_t));
	is.read(reinterpret_cast<char *>(&numActiveNodes), sizeof(int64_t));

	int64_t treeSize = fileSize - ((2*sizeof(Point3i)) + ( 2 *sizeof(int)) + maxTreeDepth + 1 + (3 * sizeof(int64_t)));
	std::cout << treeSize << std::endl;

	// Allocate memory for distanceMap & tree
	distanceMap.resize(maxTreeDepth + 1);
	tree.bits.resize(treeSize);

	is.read(reinterpret_cast<char *>(&distanceMap[0]), maxTreeDepth + 1);
	is.read(reinterpret_cast<char *>(&tree.bits[0]), treeSize);

	// Close the file
	is.close();

}

bool VolumeKdtree::pruneTreeRecursive(int64_t rootIdx) {

	int64_t reconIdx;
	bool meetsErrorThreshold = true;
	bool leftSub = true, rightSub = true;
	int rootDepth = (int)floor(log2(rootIdx + 1));

	// prune subtrees
	if (rootDepth < origTreeDepth) {
		if (parallelism) {
		//if (false) {
			parallel_invoke(
				[this, rootIdx, &leftSub] { leftSub = this->pruneTreeRecursive(2 * rootIdx + 1); },
				[this, rootIdx, &rightSub] { rightSub = this->pruneTreeRecursive(2 * rootIdx + 2); }
			);
		}
		else {
			leftSub = pruneTreeRecursive(2 * rootIdx + 1);
			rightSub = pruneTreeRecursive(2 * rootIdx + 2);
		}
	}

	if (rootDepth == origTreeDepth) {
		reconIdx = rootIdx - firstOrigLeaf;
		meetsErrorThreshold = abs(recon[reconIdx] - temp[reconIdx]) < tolerance;
	}

	// prune root, return whether or not this node was pruned
	if (leftSub && rightSub && (tree[rootIdx] == 0) && meetsErrorThreshold) {
		tree[rootIdx] = 3;
		return true;
	}
	return false;
}

void VolumeKdtree::convertToPreorder() {

	TwoBitArray preorderTree;
	preorderTree.resize(numMaxNodes);

	std::stack<std::tuple<int64_t, int, bool, long long>> stack; // node input index, node depth, eval, zeroStartIdx
	int64_t inputIdx, outputIdx, reconIdx, zeroStartIdx; // index into 'tree', index into 'preorderTree', index into 'recon'/'temp', starting idx into 'preorderTree' for pruning
	bool eval; // whether to run encoder
	int code, depth; // 2-bit node code from 'tree', node depth

	// Initialize output index & stack
	outputIdx = 0;
	stack.push(std::make_tuple(0, 0, false, -1));

	while (stack.empty() == false) {

		// Get info of node on top of stack
		inputIdx = std::get<0>(stack.top());
		depth = (int)std::get<1>(stack.top());
		eval = std::get<2>(stack.top()); 
		zeroStartIdx = std::get<3>(stack.top());
		code = (int)tree[inputIdx];

		// If current node is a leaf or a node along grown branch, check if need to run encoder
		if (depth >= origTreeDepth) {
			reconIdx = inputIdx - firstOrigLeaf;
			if (eval) {
				//recon[reconIdx] = encodeNode(inputIdx, depth, recon[reconIdx], true, temp[reconIdx]); // tree[inputIdx] is replaced with new code of this node
				recon[reconIdx] = encodeNode(inputIdx, recon[reconIdx], distanceMap[depth], true, temp[reconIdx]);
				code = (int)tree[inputIdx];

				if (zeroStartIdx != -1) { // zero present in grown branch
					if (code != 0)
						zeroStartIdx = -1;
				}
				else { // zero NOT present in grown branch
					if (code == 0)
						zeroStartIdx = outputIdx;
				}
			}
			else {
				if (depth > origTreeDepth)
					code = 3;
			}
		}
		
		// Insert current node into 'preorderTree'
		preorderTree[outputIdx++] = code;

		// Pop current node from stack
		stack.pop();

		// Continue loop if at maximum tree depth
		if (depth >= maxTreeDepth || code == 3) {
			// Prune tree if zeroStartIdx indicated
			if (zeroStartIdx != -1)
				for (int64_t i = zeroStartIdx; i < outputIdx; i++)
					preorderTree[i] = 3;
			continue;
		}

		// If at or below original tree depth, check node error
		if (depth >= origTreeDepth) {
			// If error greater than tolerance, grow branch by adding left child to stack
			if (abs(recon[reconIdx] - temp[reconIdx]) > tolerance) {
				stack.push(std::make_tuple(inputIdx, depth + 1, true, zeroStartIdx));
				continue;
			}
			// Otherwise, add '3' code & continue loop
			else {
				stack.push(std::make_tuple(inputIdx, depth + 1, false, zeroStartIdx)); 
				continue;
			}
				
		}

		// If above original tree depth, push both children to stack
		stack.push(std::make_tuple(2 * inputIdx + 2, depth + 1, false, zeroStartIdx)); // right
		stack.push(std::make_tuple(2 * inputIdx + 1, depth + 1, false, zeroStartIdx)); // left
	}

	// Swap preorderTree with tree
	// tree will size of numActiveNodes
	numActiveNodes = outputIdx;
	preorderTree.resize(numActiveNodes);
	tree.swap(preorderTree);
	preorderTree.clear();
	preorderTree.shrink_to_fit();

	//temp.clear();
	//temp.shrink_to_fit();
	//recon.clear();
	//recon.shrink_to_fit();
}

void VolumeKdtree::levelCut(int cutDepth, std::vector<byte> &outData) {

	// Set global variables
	output = &outData;
	queryDepth = cutDepth;
	  
	// Resize output vector
	output->resize(X*Y*Z);

	// Processing variables
	std::stack<std::tuple<int64_t, byte, byte, Point3i, Point3i>> stack; //idx, depth, scalar, minBound, maxBound
	int64_t idx, nextRight, nextLeft, numCells, x, y, z, c;
	int depth, code, splitDim, i;
	byte scalar;
	Point3i minBound, maxBound, extent;

	// Push root to stack to start
	stack.push(std::make_tuple(0, 0, distanceMap[0], rootMin, rootMax));

	while (stack.empty() == false) {

		// Process node at top of stack
		idx = std::get<0>(stack.top()); // index into (preorder) tree
		depth = (int)std::get<1>(stack.top());
		scalar = std::get<2>(stack.top());
		minBound = std::get<3>(stack.top());
		maxBound = std::get<4>(stack.top());
		code = (int)tree[idx];

		// If this node is a leaf....
		if (code == 3 || depth == cutDepth) {

			// Populate output vector within this node's bounding box
			for (x = minBound[0]; x < maxBound[0]; x++) {
				for (y = minBound[1]; y < maxBound[1]; y++) {
					for (z = minBound[2]; z < maxBound[2]; z++) {
						c = getCell(x, y, z);
						(*output)[c] = scalar;
					}
				}
			}
			stack.pop();

			// Push next right child to stack, if it exists
			nextRight = idx + 1;
			if (nextRight < numActiveNodes) {

				// Get parent info of next right child node from top of stack
				idx = std::get<0>(stack.top());
				depth = (int)std::get<1>(stack.top());
				scalar = std::get<2>(stack.top());
				minBound = std::get<3>(stack.top());
				maxBound = std::get<4>(stack.top());

				stack.pop();

				// Find scalar of next right node
				code = (int)tree[nextRight];
				if (code == 1)
					scalar = (byte)std::min(255.0, (double)scalar + (double)distanceMap[depth + 1]);
				else if (code == 2)
					scalar = (byte)std::max(0.0, (double)scalar - (double)distanceMap[depth + 1]);

				// Find bounding box of next right node
				extent = maxBound - minBound;
				numCells = extent.prod();
				if (numCells > 1) {
					splitDim = depth % MAX_DIM;
					i = 0;
					while (extent[splitDim] == 1) {
						splitDim = (depth + ++i) % MAX_DIM;
					}
					minBound[splitDim] = (minBound[splitDim] + maxBound[splitDim]) / 2;
				}

				stack.push(std::make_tuple(nextRight, depth + 1, scalar, minBound, maxBound));
			}
		}

		// If this node is an internal node, push next left child to stack
		else {

			if (depth >= origTreeDepth)
				stack.pop();

			nextLeft = idx + 1;

			// Find scalar of next left node
			code = (int)tree[nextLeft];
			if (code == 1)
				scalar = (byte)std::min(255.0, (double)scalar + (double)distanceMap[depth + 1]);
			else if (code == 2)
				scalar = (byte)std::max(0.0, (double)scalar - (double)distanceMap[depth + 1]);

			// Find bounding box of next left node
			extent = maxBound - minBound;
			numCells = extent.prod();
			if (numCells > 1) {
				splitDim = depth % MAX_DIM;
				i = 0;
				while (extent[splitDim] == 1) {
					splitDim = (depth + ++i) % MAX_DIM;
				}
				maxBound[splitDim] = (minBound[splitDim] + maxBound[splitDim]) / 2;
			}
			
			stack.push(std::make_tuple(nextLeft, depth + 1, scalar, minBound, maxBound));
		}
	}
}