#ifndef MID_RANGE_TREE_H
#define MID_RANGE_TREE_H

#include "TwoBitArray.h"
#include <vector>
#include <string>
#include <math.h>
#include <algorithm> 
#include <fstream>
#include <stack>
#include "point.h"
#include <numeric>
#include <random>
#include <ctime>
//#include <omp.h>
#include <ppl.h>
#include "DebugTimer.h"

using namespace concurrency;


typedef unsigned char byte;
typedef TPoint<int64_t, 3> Point3i;
//typedef byte MinMax[2];



struct MinMax {
	byte min;
	byte max;

	MinMax() {
		min = 0;
		max = 0;
	}

	MinMax(byte n, byte x) {
		min = n;
		max = x;
	}

};

/**
*	Class MidRangeTree
*	Constructs a 3D Kd-tree from input volume dataset.
*	Tree is stored in an implicit 2-bit array representing
*   progressively compressed scalar values.
*
*/
class MidRangeTree {
public:

	/************************** Public Members ****************************/

	/* Root Bounding Box */
	Point3i rootMin;
	Point3i rootMax;

	/* 2-bit nodes represent change in scalar value from parent to child */
	TwoBitArray tree;
	TwoBitArray tree_range;

	/* Map from tree depth to distance in scalar values from parent to child */
	std::vector<byte> distanceMap;
	std::vector<byte> distanceMap_range;

	/* Depth of tree (does not include root level) */
	int maxTreeDepth;
	int origTreeDepth;

	/* Number of non-null nodes in tree */
	int64_t numActiveNodes;

	/* Dimensions of original volume data */
	int64_t X, Y, Z;

	/* Knobs to adjust amount of extra levels to add to
	tree for increased accuracy */
	int tolerance;
	int maxEpochs;

	/* Buffers for queries */
	std::vector<byte> * output;
	std::vector<byte> * output2;
	int queryDepth;

	/*********************** Public Functions ***************************/

	/* Default Constructor */
	MidRangeTree() {
		rootMin = 0;
		rootMax = 0;
		tolerance = 6;
		maxEpochs = 5;
	}

	/** Constructor

	@param inData - pointer to original volume dataset as 1-D array
	@param x - number of cells along x-axis of dataset
	@param y - number of cells along y-axis of dataset
	@param z - number of cells along z-axis of dataset
	*/
	MidRangeTree(std::vector<byte> &inData, int64_t x, int64_t y, int64_t z) {
		data = &inData;
		X = x;
		Y = y;
		Z = z;
		rootMin = 0;
		rootMax = { X,Y,Z };
		tolerance = 6;
		maxEpochs = 5;
	}

	/* Deconstructor */
	~MidRangeTree() {
		std::vector<byte>().swap(distanceMap);
		std::vector<byte>().swap(temp);
		tree.bits.clear();
		tree.bits.shrink_to_fit();
		tree_range.bits.clear();
		tree_range.bits.shrink_to_fit();
	}

	/**
	Builds Kd-tree from input volume dataset.
	*/
	void build(bool useThreads = true);

	/**
	Returns data of nodes at specific level of the tree.
	This function traverses the depth-first tree array.

	@param cutDepth - desired level of tree.
	@param outData - pointer to output buffer.
	*/
	void levelCut(int cutDepth, std::vector<byte> &outData);

	/**
	Calculate maximum cell error.

	@return maximum L1 error between query data and original dataset.
	*/
	int measureMaxError();

	/**
	Calculates mean cell error.

	@return mean L1 error between query data and original dataset.
	*/
	double measureMeanError();

	void queryError(std::vector<byte> &outData);

	/**
	Writes class to binary file.

	@param filename Path of output binary file.
	*/
	void save(std::string filename);

	/**
	Opens serialized class from binary file

	@param filename Path of input binary file.
	*/
	void open(std::string filename);

	/**
	Sets the error tolerance for nodes.
	Default is 6.

	@param Error tolerance.
	*/
	void setErrorTolerance(int errorTolerance);

	/**
	Sets the maximum amount of epochs when finding distance map values
	using gradient descent.
	Default is 8.

	@param Number of epochs.
	*/
	void setMaxEpochs(int epochs);

	void convertToByteArray(std::vector<byte> &byteArray);


private:

	/************************* Private Members ********************************/

	/* Temporary tree for construction (before compression) */
	std::vector<byte> temp;
	std::vector<byte> temp_range;

	/* Temporary vector of reconstructed nodes */
	std::vector<byte> recon;
	std::vector<byte> recon_range;

	/* Pointer to original volume data */
	std::vector<byte> * data;

	/* 3-D tree */
	static int const MAX_DIM = 3;

	/* Temporary variables for adding levels to tree */
	int64_t numOrigNodes;
	int64_t numMaxNodes;
	int64_t firstOrigLeaf;

	bool parallelism;


	/*********************** Private Functions *******************************/

	/**
	Recursive tree building function.
	Populates a true byte node and recurses on node's children.

	@param idx Node index.
	@param depth Depth of node in the tree.
	@param minBound 3D point defining minimum coordinate of node's bounding box.
	@param maxBound 3D point defining maximum coordinate of node's bounding box.

	@return Node's bounding box.
	*/
	MinMax buildRecursive(int64_t idx, int depth, Point3i minBound, Point3i maxBound);

	/**
	Compresses 8-bit trees to 2-bit trees using gradient descent to 
	find optimal distance map values.
	*/
	void compressGradientDescent();

	void compressGradientDescent_range();

	/**
	Return cell index into 1-D buffer of 3-D volume.

	@param x x-coordinate.
	@param y y-coordinate.
	@param z z-coordinate.

	@return cell index into data array
	*/
	inline int64_t getCell(int64_t x, int64_t y, int64_t z);

	/**
	Recursively prunes subtrees by setting node code to 3

	@param rootIdx Root index of subtree to prune (index into breath-first array).
	*/
	bool pruneTreeRecursive(int64_t rootIdx);

	/**
	Converts breath-first full tree array into depth-first unbalanced array representation.
	*/
	void convertToPreorder();

	byte encodeNodeEstimate(int64_t idx, byte compressedParent, double * estimateSum, double * estimateCount);

	byte encodeNodeEstimate_range(int64_t idx, byte compressedParent, double * estimateSum, double * estimateCount);

	byte encodeNode(int64_t idx, byte compressedParent, byte distanceVal, bool fillTree = true,
		double nodeTruth = -1.0, double * returnError = nullptr);

	byte encodeNode_range(int64_t idx, byte compressedParent, byte distanceVal, bool fillTree = true,
		double nodeTruth = -1.0, double * returnError = nullptr);
};

#endif MID_RANGE_TREE_H
