#ifndef UNITBRICK_H
#define UNITBRICK_H

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLM Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

/**
* Class UnitBrick
* Stores a unit brick (1x1x1/2) on GPU surrounding origin using a vertex buffer object.
**/
class UnitBrick {

public:

	/**** Public Attributes *****/

	/* Buffer ID's */
	GLuint cubeVBO, cubeVAO, cubeIndsVBO;



	/***** Public Functions *****/

	/* Default constructor */
	UnitBrick() {}

	/* Default deconstructor */
	~UnitBrick() {}

	/* Creates and stores the unit cube on the GPU using VBOs*/
	void Setup() {

		// Generate unit cube vertex array and vertex buffer objects
		glGenVertexArrays(1, &cubeVAO);
		glGenBuffers(1, &cubeVBO);
		glGenBuffers(1, &cubeIndsVBO);

		// Unit cube vertices 
		glm::vec3 vertices[8] = { glm::vec3(-0.5f,-0.5f,-0.25f),
			glm::vec3(0.5f,-0.5f,-0.25f),
			glm::vec3(0.5f, 0.5f,-0.25f),
			glm::vec3(-0.5f, 0.5f,-0.25f),
			glm::vec3(-0.5f,-0.5f, 0.25f),
			glm::vec3(0.5f,-0.5f, 0.25f),
			glm::vec3(0.5f, 0.5f, 0.25f),
			glm::vec3(-0.5f, 0.5f, 0.25f) };

		// Unit cube indices
		GLushort cubeIndices[36] = { 0,5,4,
			5,0,1,
			3,7,6,
			3,6,2,
			7,4,6,
			6,4,5,
			2,1,3,
			3,1,0,
			3,0,7,
			7,0,4,
			6,5,2,
			2,5,1 };

		// Bind cube array and buffer
		glBindVertexArray(cubeVAO);
		glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);

		// Pass cube vertices to buffer object memory
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), &(vertices[0].x), GL_STATIC_DRAW);

		// Enable vertex attribute array for position
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		// Pass indices to element array  buffer
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeIndsVBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cubeIndices), &cubeIndices[0], GL_STATIC_DRAW);

		// Unbind array
		//glBindVertexArray(0);
		Unbind();
	}

	/* Method to draw the unit cube */
	void Draw() {
		glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);
	}

	void Bind() {
		glBindVertexArray(cubeVAO);
		glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeIndsVBO);
	}

	void Unbind() {
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	void Delete() {
		glDeleteVertexArrays(1, &cubeVAO);
		glDeleteBuffers(1, &cubeVBO);
		glDeleteBuffers(1, &cubeIndsVBO);
	}
};
#endif //UNITBRICK_H