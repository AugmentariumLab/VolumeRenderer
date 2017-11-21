#ifndef OCTREE_H
#define OCTREE_H

#include "TwoBitArray.h"
#include <vector>
#include <string>
#include <math.h>
#include <algorithm> 


typedef unsigned char byte;
typedef float point3D[3];
/*
struct IdxRange {
	int64_t first;
	int64_t last;

	IdxRange() : first(0), last(0) { }

	IdxRange(int64_t f, int64_t l) : first(f), last(l) {}
};
*/

struct BoundingBox {
	point3D min;
	point3D max;

	BoundingBox() {
		for (int i = 0; i < 3; i++) {
			min[i] = 0;
			max[i] = 0;
		}
	}

	// Create a collapsed bounding box from a single point
	BoundingBox(const point3D &p) {
		memcpy(min, p, sizeof(point3D));
		memcpy(max, p, sizeof(point3D));
	}

	// Create a bounding box from two points
	BoundingBox(const point3D &minP, const point3D &maxP) {
		memcpy(min, minP, sizeof(point3D));
		memcpy(max, maxP, sizeof(point3D));
	}

	std::vector<BoundingBox> split2(int axis) {
		float mid = (min[axis] + max[axis]) / 2.0f;

		point3D _max, _min;
		memcpy(_max, max, sizeof(point3D));
		_max[axis] = mid;
		memcpy(_min, min, sizeof(point3D));
		_min[axis] = mid;

		std::vector<BoundingBox> children (2);
		children[0] = BoundingBox(min, _max);
		children[1] = BoundingBox(_min, max);

		return children;
	}

	bool isCube() {
		double extent = max[0] - min[0];
		if (max[1] - min[1] != extent)
			return false;
		if (max[2] - min[2] != extent)
			return false;

		return true;
	}
	// returns extent if is an even cube, -1 if not
	void getExtent(point3D extents) {
		for (int i = 0; i < 3; i++) {
			extents[i] = max[i] - min[i];
		}
	}

	// Return the position of a bounding box corner
	void getCorner(int index, point3D &corner) {
		for (int i = 0; i<3; ++i)
			corner[i] = (index & (1 << i)) ? max[i] : min[i];
	}

	int64_t getVolume() {
		point3D e;
		getExtent(e);
		return e[0] * e[1] * e[2];
	}

	std::vector<BoundingBox> split8() {

		// return empty vector if box is not an even cube
		if (!isCube())
			return std::vector<BoundingBox>();

		float childExtent = (max[0] - min[0]) / 2.0f;

		int i;
		point3D mid; 
		for (i = 0; i < 3; i++) {
			mid[i] = (min[i] + max[i]) / 2.0f;
		}

		// loop thru corners of mid -> min box 
		BoundingBox firstChild = BoundingBox(min, mid);
		std::vector<BoundingBox> children(8);
		children[0] = firstChild;

		for (i = 1; i < 8; i++) {
			point3D _min, _max;
			firstChild.getCorner(i, _min);
			for (int j = 0; j < 3; j++)
				_max[j] = _min[j] + childExtent;

			children[i] = BoundingBox(_min, _max);
		}
		
		return children;
	}

	// Print a string representation of the bounding box
	void print() const {
		std::cout << "Box[min= (" << min[0] << " " << min[1] << " " << min[2] << "), max= (";
		std::cout << max[0] << " " << max[1] << " " << max[2] << ") ]" << std::endl;
	}
};

#pragma pack(push, 1)
struct TempNode {
	byte parentScalar;
	byte scalar;

	TempNode() : parentScalar(0), scalar(0) {}

	TempNode(byte pS, byte s) : parentScalar(pS), scalar(s) {}
};
#pragma pack(pop)

/**
*	Class Octree
*	Constructs an octree from input volume dataset.
*	Octree is uses implicit 2-bit array. If volume is not uniform in all directions,
*	it will start as a binary tree then transform to an octree.
*
*/
class Octree {
public:

	


	/***** Public Members *****/

	/* Root Bounding Box */
	BoundingBox rootBox;

	/* 2-bit nodes represent change in scalar value from parent to child */
	TwoBitArray * tree;

	/* Depth in the tree where nodes start to have 8 children. 
	Nodes at depths under this number have 2 children*/
	int octreeStartDepth;
	
	/* Map from tree depth to change in scalar value */
	std::vector<byte> deltaScalarMap;

	/* Which axes are split in binary part of tree*/
	std::vector<int> splitDims;




	/***** Public Functions *****/

	/* Constructor */
	Octree() {;}

	/* Deconstructor */
	~Octree() {;}

	/* Builds octree from input volume dataset */
	void build(std::vector<byte> * inData, int64_t X, int64_t Y, int64_t Z);

	/* Writes class to binary file */
	//void save(std::string filename);

	/* Opens serialized class from binary file */
	//void open(std::string filename); 


private:

	/***** Private Members *****/

	std::vector<TempNode> tempTree;

	/* Pointer to original volume data */
	std::vector<byte> *data;

	/* 3-D vector version of original volume data */
	std::vector<std::vector<std::vector<byte>>> data3D;

	/* Dimensions of original volume data */
	int64_t X, Y, Z;



	/***** Private Functions *****/

	/* Recursive tree building function, called by build(data) */
	void buildRecursive(int64_t idx, int depth, BoundingBox box, byte parentScalar);

	//byte encodeDeltaScalar(int64_t idx, int depth, byte parentScalar, byte outScalar, bool useMap);

	static int calcOctreeDepth(int64_t numCells);

	inline int64_t getCell(int64_t x, int64_t y, int64_t z);

	static int64_t calcNumOctreeNodes(int depth);

};

#endif OCTREE_H
