#ifndef HASHED_KDTREE_H
#define HASHED_KDTREE_H

#include <vector>
#include <string>
#include <math.h>
#include <algorithm> 
#include <fstream>
#include "point.h"
#include "TwoBitArray.h"
#include <random>
#include <unordered_map>



typedef unsigned char byte;
typedef TPoint<int64_t, 3> Point3i;

/**
*	Class HashedKdtree
*	Constructs an unbalanced 3D kd-tree from an input volume dataset.
*	Tree is stored in a hash table where the key is a morton code
*	and the values represent progressively compressed scalar values.
*
*/
class HashedKdtree {
public:

	/***** Public Members *****/

	/* Root Bounding Box */
	Point3i rootMin;
	Point3i rootMax;

	/* 2-bit nodes represent change in scalar value from parent to child */
	TwoBitArray treeData; //0=nochange; 1=poschange; 2=negchange; 3=collision
	TwoBitArray treeStructure; //0=no children; 1=left child only; 2=right child only; 3=both children
	std::unordered_map<int64_t, int32_t> collisionNodes; // morton code -> idx to collision vectors
	TwoBitArray treeDataCollisions; // tree data overflow vector for morton code collisions
	TwoBitArray treeStructureCollisions; // tree structure overflow vector for morton code collisions

	/* Map from tree depth to distance in scalar values from parent to child */
	std::vector<byte> distanceMap;

	/* Depth of tree (does not include root level) */
	int treeDepth;

	int64_t numCollisions;
	
	int64_t hashMask;

	/* Dimensions of original volume data */
	int64_t X, Y, Z;

	/* Buffers for queries */
	std::vector<byte> * output;
	std::vector<byte> * output2;
	int queryDepth;

	/* Leaf accuracy tolerance */
	int tolerance;


	/***** Public Functions *****/

	/* Default Constructor */
	HashedKdtree() {
		rootMin = 0;
		rootMax = 0;
	}

	HashedKdtree(std::vector<byte> &inData, int64_t x, int64_t y, int64_t z) {
		data = &inData;
		X = x;
		Y = y;
		Z = z;
		rootMin = 0;
		rootMax = { X,Y,Z };
		tolerance = 4;
		addLevelDistance = 64;
		maxAddLevels = 4;
		std::random_device random_dev;
		std::mt19937 generator(random_dev());
		gen = generator;
		lastCollisionIdx = 0;
	}

	/* Deconstructor */
	~HashedKdtree() {
		std::vector<double>().swap(distanceSums);
		std::vector<double>().swap(distanceCounts);
		std::vector<byte>().swap(distanceMap);
		std::vector<byte>().swap(temp);
		//delete &tree;
	}

	/**
	Builds octree from input volume dataset.

	@param inData - pointer to original volume dataset in 1-D array
	@param X - number of cells along x-axis of dataset
	@param Y - number of cells along y-axis of dataset
	@param Z - number of cells along z-axis of dataset
	*/
	void build();

	/**
	Return nodes at specific level of the tree.

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


private:

	/***** Private Members *****/

	int origTreeDepth;
	int64_t numNodes;

	/* Temporary arrays for contruction (before compression) */
	std::vector<byte> temp; // true node scalar value
	std::vector<byte> tempCollisions; // true node scalar value overflow for collisions
	std::vector<int64_t> visited; // temporary array of morton codes (hashed key -> morton code)
	

	/* Pointer to original volume data */
	std::vector<byte> * data;

	/* 3-D tree */
	static int const MAX_DIM = 3;

	/* Temporary vectors for distanceMap construction */
	std::vector<double> distanceSums;
	std::vector<double> distanceCounts;

	/* */
	byte addLevelDistance;
	int maxAddLevels;
	std::vector<int> childVec = { 1,2 };
	std::mt19937 gen;
	int64_t lastCollisionIdx;


	/***** Private Functions *****/

	/**
	Recursive tree building function.
	Populates a node and recurses on node's children.

	@param idx Node index.
	@param depth Depth of node in the tree.
	@param minBound 3D point defining minimum coordinate of node's bounding box.
	@param maxBound 3D point defining maximum coordinate of node's bounding box.
	@param decodedParent scalar value of node's parent as estimated by decoder.
	*/
	void buildRecursive(int64_t mcode, int depth, Point3i minBound, Point3i maxBound, byte decodedParent);

	void levelCutRecursive(int64_t mcode, int depth, Point3i minBound, Point3i maxBound, byte parentEstimate);

	/**
	Assigns 2-bit code to tree node & returns decoded scalar value for that node.

	@param idx Node Index.
	@param depth Depth of node in the tree.
	@param decodedParent scalar value of node's parent as estimated by decoder.
	@param useMap True is using populated distanceMap.

	@return Node's scalar value as estimated by decoder.
	*/
	byte encodeNode(int64_t mcode, int depth, byte decodedParent, byte trueNode, bool useMap, int collisionIdx=-1);

	/**
	Return cell index into 1-D buffer of 3-D volume.

	@param x x-coordinate.
	@param y y-coordinate.
	@param z z-coordinate.

	@return cell index into data array
	*/
	inline int64_t getCell(int64_t x, int64_t y, int64_t z);

	inline int64_t hash(int64_t code);


	/**
	Recursive compression function.
	Assigns a delta scalar code value to a node.

	@param idx Node Index.
	@param depth Depth of node in the tree.
	*/
	void compressTreeRecursive(int64_t mcode, int depth, byte compressedParent, int trueNodeOverride=-1);
};

#endif VOLUME_OCTREE_H
