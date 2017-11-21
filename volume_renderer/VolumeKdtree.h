#ifndef VOLUME_KDTREE_H
#define VOLUME_KDTREE_H

#include "TwoBitArray.h"
#include <vector>
#include <string>
#include <math.h>
#include <algorithm> 
#include <fstream>
#include "point.h"


typedef unsigned char byte;
typedef TPoint<float, 3> Point3f;
typedef TPoint<int64_t, 3> Point3i;

/**
	Struct BoundingBox.
	A bounding box for cells in volume dataset.
*/
struct BoundingBox {
	Point3f min; // minimum point defining bounding box
	Point3f max; // maximum point defining bounding box
	Point3i volMin; // minimum coordinate of volume 
	Point3i volMax; // maximum coordinate of volume 

	// Construct a collapsed default box at (0,0,0)
	BoundingBox() : min(0.0f), max(0.0f), volMin(0), volMax(0) { }

	// Construct a bounding box from two points
	BoundingBox(const Point3f &minP, const Point3f &maxP, 
		const Point3i &volMin, const Point3i &volMax) : 
		min(minP), max(maxP), volMin(volMin), volMax(volMax) { }

	// Returns X,Y,Z extents of bounding box
	Point3f getExtents() {
		return max - min;
	}

	// Return the position of a bounding box corner
	Point3f getCorner(int index) {
		Point3f corner;
		for (int i = 0; i<3; ++i)
			corner[i] = (index & (1 << i)) ? max[i] : min[i];
		return corner;
	}

	// Return the volume of the bounding box
	float getVolume() {
		return (max - min).prod();
	}

	// Return number of cells the box encompasses in original volume
	int64_t getNumCells() {
		return (volMax- volMin + Point3i(1, 1, 1)).prod();
	}

	// Return the center point of the bounding box 
	Point3f getCenter() {
		return (max + min) * 0.5f;
	}

	// Print a string representation of the bounding box
	void print() {
		std::cout << "Box[min= (" << min[0] << " " << min[1] << " " << min[2] << "), max= (";
		std::cout << max[0] << " " << max[1] << " " << max[2] << ")" << std::endl;
		std::cout << "volMin = (" << volMin[0] << " " << volMin[1] << " " << volMin[2] << "), volMax = (";
		std::cout << volMax[0] << " " << volMax[1] << " " << volMax[2] << ")" << std::endl;
	}
};

#pragma pack(push, 1)
struct TempNode {
	byte parentScalar; // scalar value of node's parent
	byte scalar; // average scalar value in bounding box

	TempNode() : parentScalar(0), scalar(0) {}

	TempNode(byte pS, byte s) : parentScalar(pS), scalar(s) {}
};
#pragma pack(pop)


/**
*	Class VolumeKdtree
*	Constructs a 3D kd-tree from input volume dataset.
*	Tree is stored in an implicit 2-bit array representing compressed scalar values. 
*
*/
class VolumeKdtree {
public:

	/***** Public Members *****/

	/* Root Bounding Box */
	BoundingBox rootBox;

	/* 2-bit nodes represent change in scalar value from parent to child */
	TwoBitArray tree;

	/* Map from tree depth to distance in scalar values from parent to child */
	std::vector<byte> distanceMap;

	/* Depth of tree (does not count root) */
	int treeDepth;

	/* Depth of tree in which a node represents 1 cell in volume */
	//int oneCellDepth;

	/* Dimensions of original volume data */
	int64_t X, Y, Z;


	double mse;
	double maxError;

	/***** Public Functions *****/

	/* Constructor */
	VolumeKdtree() { 
		added = false;}

	/* Deconstructor */
	~VolumeKdtree() {
		std::vector<byte>().swap(distanceMap);
		std::vector<double>().swap(distanceSums);
		std::vector<double>().swap(nonzeroNodeCounts);
		std::vector<TempNode>().swap(tempTree);
		//free(tree);
		tree = NULL;
	}

	/**
		Builds octree from input volume dataset.

		@param inData - pointer to original volume dataset in 1-D array
		@param X - number of cells along x-axis of dataset
		@param Y - number of cells along y-axis of dataset
		@param Z - number of cells along z-axis of dataset
	*/
	void build(std::vector<byte> &inData, int64_t X, int64_t Y, int64_t Z);

	void levelCut(int cutDepth, std::vector<byte> &outData);

	/* Writes class to binary file */
	void save(std::string filename);

	/* Opens serialized class from binary file */
	void open(std::string filename); 


private:

	/***** Private Members *****/

	bool added;

	/* Temporary tree for contruction (before compression) */
	std::vector<TempNode> tempTree;

	/* Pointer to original volume data */
	std::vector<byte> * data;

	/* Map from depth to sum of distances from parents */
	std::vector<double> distanceSums;

	/* Map from depth to number of non-zero nodes 
	(nodes that have some significant change from parent) */
	std::vector<double> nonzeroNodeCounts;

	/* 3-D tree */
	static int const MAX_DIM = 3;


	int64_t numZeroNodes;

	/***** Private Functions *****/

	/** 
		Recursive tree building function, called by build() function.
		Populates a node and recurses on node's children.
		
		@param idx - Node index
		@param depth - Node depth in tree
		@param box - Bounding box of node
		@param parentScalar - scalar value of Node's parent
	*/
	void buildRecursive(int64_t idx, int depth, BoundingBox box, byte parentScalar);

	/**
		Encodes 2-bit delta scalar value. 
		Delta Codes:
			0 = no change
			1 = add distance
			2 = subtract distance
			3 = no change & no children

		@param idx Node index.
		@param depth Node depth in tree.
		@param parent Scalar value of Node's parent.
		@param estimate Estimate of Node's scalar value, as encoded by 2-bit delta (out value).
		@param useMap True if using a populated distanceMap.

		@return 2-bit delta code
	*/
	byte encodeDeltaScalar(int64_t idx, int depth, byte parent, byte &estimate, bool useMap);

	/**
		Return cell index into 1-D buffer of 3-D volume.

		@param x x-coordinate.
		@param y y-coordinate.
		@param z z-coordinate.

		@return cell index into data array
	*/
	inline int64_t getCell(int64_t x, int64_t y, int64_t z);

	/**
		Convert tree from temporary 2-byte representation to 2-bit representation.
	*/
	void compressTree(bool computeMap);

	/** 
		Recursive compression function, called from compressTree().
		Assigns a delta scalar code value to a node.

		@param idx Node Index.
		@param depth Depth of node in the tree.
	*/
	void compressTreeRecursive(int64_t idx, int depth);

	void levelCutRecursive(int64_t idx, int depth, byte parentScalar, BoundingBox box, 
		int cutDepth, std::vector<double> &sumData, std::vector<double> &countData);

	//void addLevels(int numLevels);

	void buildFromLeaves(int64_t idx, int depth, byte parentScalar, BoundingBox box, int leafDepth);

	/* Prunes compressed tree using the 3 code. */
	void pruneTree();
};

#endif VOLUME_OCTREE_H
