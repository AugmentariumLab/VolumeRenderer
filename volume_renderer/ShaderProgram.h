#ifndef SHADERPROGRAM_H
#define SHADERPROGRAM_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

class ShaderProgram
{
public:
	/****** Public Attributes *****/

	/* ID's to shaders */
	GLuint shaderIds[3]; //0->vertexshader, 1->fragmentshader, 2->geometryshader

	/* ID to shader program */
	GLuint programId;

	/* Uniform variables */
	std::map<std::string, GLuint> uniforms;

	/* Attribute variables */
	std::map<std::string, GLuint> attributes;


	/****** Public Functions ******/

	/* Default Constructor */
	ShaderProgram() {
		programId = 0;
		for (int i = 0; i < 3; i++)
			shaderIds[i] = 0;
	}

	/* Default Deconstructor */
	~ShaderProgram() {
		uniforms.clear();
	}

	/* Loads a shader from a file */
	bool loadShaderFromFile(const char* path, GLenum shaderType)
	{
		//Open file
		GLuint shaderID = 0;
		std::string shaderString;
		std::ifstream sourceFile(path);

		//Source file loaded
		if (sourceFile)
		{
			//Get shader source
			shaderString.assign((std::istreambuf_iterator< char >(sourceFile)), std::istreambuf_iterator< char >());

			//Create shader ID
			shaderID = glCreateShader(shaderType);

			//Set shader source
			const GLchar* shaderSource = shaderString.c_str();
			glShaderSource(shaderID, 1, (const GLchar**)&shaderSource, NULL);

			//Compile shader source
			glCompileShader(shaderID);

			//Check shader for errors
			GLint shaderCompiled = GL_FALSE;
			glGetShaderiv(shaderID, GL_COMPILE_STATUS, &shaderCompiled);
			if (shaderCompiled != GL_TRUE)
			{
				printf("Unable to compile shader %d!\n\nSource:\n%s\n", shaderID, shaderSource);
				glDeleteShader(shaderID);
				shaderID = 0;
				return false;
			}

			// Save shader ID
			if (shaderType == GL_VERTEX_SHADER)
				shaderIds[0] = shaderID;
			else if (shaderType == GL_FRAGMENT_SHADER)
				shaderIds[1] = shaderID;
			else
				shaderIds[2] = shaderID;

			return true;
		}
		else
		{
			printf("Unable to open file %s\n", path);
			return false;
		}
	}

	bool createProgram() {
		int i;

		programId = glCreateProgram();

		// Attach compiled shaders to shader program object
		for (i = 0; i < 3; i++) {
			if (shaderIds[i] != 0)
				glAttachShader(programId, shaderIds[i]);
		}

		glLinkProgram(programId);

		// Check if shader program failed
		GLint success;
		GLchar infoLog[512];
		glGetProgramiv(programId, GL_LINK_STATUS, &success);
		if (!success) {
			glGetProgramInfoLog(programId, 512, NULL, infoLog);
			std::cout << "ERROR::SHADER::PROGRAM::COMPILATION_FAILED\n" << infoLog << std::endl;
		}

		// Delete shader objects once they've been linked to program object
		for (i = 0; i < 3; i++)
			glDeleteShader(shaderIds[i]);

		return (bool)success;
	}

	void addUniform(const std::string& uniform) {
		if (programId == 0) {
			throw std::runtime_error("\nError! Cannot get uniform location, shader program does not exist!\n");
		}
		else {
			uniforms[uniform] = glGetUniformLocation(programId, uniform.c_str());
		}
	}

	void addAttribute(const std::string& attribute) {
		if (programId == 0) {
			throw std::runtime_error("\nError! Cannot get uniform location, shader program does not exist!\n");
		}
		else {
			attributes[attribute] = glGetAttribLocation(programId, attribute.c_str());
		}
	}


	void Use() {
		glUseProgram(programId);
	}

	void StopUse() {
		glUseProgram(0);
	}
};

#endif //SHADERPROGRAM_H
