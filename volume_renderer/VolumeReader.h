#ifndef VOLUMEREADER_H
#define VOLUMEREADER_H

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <stdexcept>
#include <map>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

typedef unsigned int dim3D[3]; // (X,Y,Z)


/**
* Class VolumeReader
* Loads a volume dataset broken up into multiple binary files 
* into a single 3-D Texture for volume rendering.
*
**/
template<typename T>
class VolumeReader {

public:
	
	/***** Public Attributes *****/

	/* Expected dimensions of each brick of volume dataset */
	dim3D brickDims;

	/* Function to find dataset based on brick number (arg1) & timestep (arg2) */
	std::function<std::string(int, int)> findSourceFile; 
	
	/* Texture ID */
	GLuint textureId;
	
	/* Map from brick # to I,J,K coords of brick grid*/
	std::map<int, dim3D> * brickMap;

	/***** Public Functions *****/

	/* Default Constructor */
	VolumeReader() {
		// Default class attributes
		brickDims[0] = brickDims[1] = brickDims[1] = 0;
		findSourceFile = [](int brick, int time)
		{ throw std::runtime_error("\n\nERROR! findSourceFile function not defined.\n");
		return "";  };
	}

	/* Constructor */
	VolumeReader(dim3D brick, dim3D grid, 
		std::function<std::string(int, int)> findFileFunct,
		std::map<int,dim3D> * bMap) {

		// Define class attributes from input
		memcpy(brickDims, brick, sizeof(dim3D));
		findSourceFile = findFileFunct;
		brickMap = bMap;
	}

	/* Deconstructor */
	~VolumeReader() {
		std::vector<T>().swap(tempBrick);
		std::vector<T>().swap(data);
	}

	/**
		Load single brick from binary file to 3D texture.

		@param brick the brick number to load.
		@param timestep the timestep to load.
		@return true if brick was successfully loaded
	*/
	bool LoadBrickToTexture(int brick, int timestep) {
		
		// load brick data to 'tempBrick' vector
		bool loadSuccess = LoadVolumeFromBinaryFile(findSourceFile(brick, timestep));

		if (loadSuccess) {
			// Create OpenGL 3D Texture
			glGenTextures(1, &textureId);
			glBindTexture(GL_TEXTURE_3D, textureId);

			// Set the texture parameters
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

			// Load data vector into 3D Texture
			glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, (GLsizei)brickDims[0], (GLsizei)brickDims[1], (GLsizei)brickDims[2], 0, GL_RED, GL_UNSIGNED_BYTE, &tempBrick[0]);

			// Generate Mipmaps
			glGenerateMipmap(GL_TEXTURE_3D);

			// Deallocate temporary vector
			std::vector<T>().swap(tempBrick);
		}
		else {
			std::cout << "ERROR! Texture load failure!" << std::endl;
		}

		return loadSuccess;
	}

	/** 
		Load multiple bricks from multiple binary files to single 3D texture.

		@param numBricks the number of bricks to load.
		@param I the number of bricks in the X direction of the brick grid.
		@param J the number of bricks in the Y direction of the brick grid.
		@param K the number of bricks in the Z direction of the brick grid.
		@param timestep the timetep to load.
		@return true if all bricks were successfully loaded.

	*/
	bool LoadBricksToTexture(int64_t numBricks, int64_t I, int64_t J, int64_t K, int timestep) {

		// For easy reference
		int64_t X = brickDims[0];
		int64_t Y = brickDims[1];
		int64_t Z = brickDims[2];
		int64_t XY = X*Y;
		int64_t XYZ = XY * Z;
		int64_t XYZIJ = XYZ * I * J;
		int64_t XYI = XY * I;
		int64_t XI = X * I;
		int64_t XYIJ = XY * I * J;

		// Make sure every file load is successful
		bool success = true;

		// Resize volume vector
		data.resize(XYZ * numBricks, 0);
	
		// Loop through bricks & load their data
		unsigned int globalBrickShift, globalZShift, brickZShift, globalStartIdx, brickStartIdx, brickEndIdx;
		for (int b = 0; b < numBricks; b++){

			// For easy reference, (i,j,k) coordinate of the brick
			unsigned int i = (*brickMap)[b][0];
			unsigned int j = (*brickMap)[b][1];
			unsigned int k = (*brickMap)[b][2];

			// Load brick data
			success *= LoadVolumeFromBinaryFile(findSourceFile(b, timestep));

			if (success) {
				// Calculate shift in global 'data' array index from (i,j,k) coordinate
				globalBrickShift = (k * XYZIJ) + (j * XYI) + (i * X);

				// Loop through box rows
				for (unsigned int z = 0; z < brickDims[2]; z++) {

					// Calculate shift in global 'data' array index from (z) position
					globalZShift = z * XYIJ;
					// Calculate shift in tempBrick array index from (z) position
					brickZShift = z * XY;

					for (unsigned int y = 0; y < brickDims[1]; y++) {

						globalStartIdx = globalBrickShift + globalZShift + (y * XI); //idx of 'data' to place start of row
						brickStartIdx = brickZShift + (y * X); // idx from tempBrick of the start of row
						brickEndIdx = brickStartIdx + X; // idx of tempBrick of the end of row

						// Copy row from temporary brick vector to global 'data' vector in appropriate vector range
						std::copy(tempBrick.begin() + brickStartIdx,
							tempBrick.begin() + brickEndIdx,
							data.begin() + globalStartIdx);
					}
				}
			}
			else {
				std::cout << "Load error. Brick loading terminated." << std::endl;
				return success;
			}		
		}

		// Find new volume dimensions
		GLuint64 volumeX = I * X;
		GLuint64 volumeY = J * Y;
		GLuint64 volumeZ = K * Z;
		std::cout << "TEXTURE SIZE: " << (double)(volumeX * volumeY * volumeZ) / 1e9 << " GB" << std::endl;

		// Create OpenGL 3D Texture
		glGenTextures(1, &textureId);
		glBindTexture(GL_TEXTURE_3D, textureId);

		// Set the texture parameters
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

		// Load data vector into 3D Texture
		glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, volumeX, volumeY, volumeZ, 0, GL_RED, GL_UNSIGNED_BYTE, &data[0]);

		// Generate Mipmaps
		//glGenerateMipmap(GL_TEXTURE_3D);

		// Deallocate temporary vectors
		std::vector<T>().swap(data);
		std::vector<T>().swap(tempBrick);

		return success;
	}

private:

	/***** Private Attributes *****/

	/* Temporary buffers */
	std::vector<T> data; // holds full dataset for transfer to GPU
	std::vector<T> tempBrick; // holds single brick

	/***** Private Functions *****/

	/** 
		Reads binary file of brick data and loads into 'tempBrick' vector.

		@param filename the full path to binary file to load.
		@return true if binary file was loaded successfully.
		
	*/
	bool LoadVolumeFromBinaryFile(std::string filename) {

		// Open the stream
		std::ifstream is(filename, std::ios::in | std::ios::binary);

		if (is.good()) {

			std::cout << "Loading data from: " << filename << "..." << std::endl;

			// Confirm file size matches volume type and dimensions
			int64_t expectedSize = sizeof(T) * brickDims[0] * brickDims[1] * brickDims[2];
			is.seekg(0, is.end);
			int64_t fileSize = is.tellg();
			is.seekg(0, is.beg);
			if (expectedSize != fileSize) {
				throw std::runtime_error("File size does not match expected dataset size!");
				return false;
			}

			// Resize vector based on volume dimensions;
			tempBrick.resize(brickDims[0] * brickDims[1] * brickDims[2]);

			// Load the data
			is.read((char*)&tempBrick[0], fileSize);

			// Check success of the load
			bool loadSuccess = false;
			if (is) {
				loadSuccess = true;
				std::cout << "\tLoad successful!" << std::endl;
			}
			else {
				std::cout << "\tLoad Error! Only " << is.gcount() << " could be read." << std::endl;
			}

			// Close the file
			is.close();

			// Return whether load was successful
			return loadSuccess;
		}
		else {
			// Could not open file - return failure status
			return false;
		}
	}
};

#endif //VOLUMEREADER_H