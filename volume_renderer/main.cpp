#include <vector>
#include <map>
#include <iterator>
#include <numeric>

#include "VolumeReader.h"
#include "ShaderProgram.h"
#include "UnitBrick.h"
#include "DebugTimer.h"
//#include "HashedKdTree.h"
#include "VolumeKdTree_recover.h"
//#include "MidRangeTree.h"

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// GLM Mathematics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Window dimensions
const GLuint WIDTH = 1600, HEIGHT = 1200;

// Starting positions
//glm::vec3 cameraPosStart = glm::vec3(0.0f, -1.5f, -1.5f);
//glm::vec3 cameraFrontStart = glm::vec3(0.0f, 1.0f, 1.0f);
//glm::vec3 cameraUpStart = glm::vec3(0.0f, 1.0f, 0.0f);
glm::vec3 cameraPosStart = glm::vec3(0.0f, 0.0f, -0.75f);
glm::vec3 cameraFrontStart = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 cameraUpStart = glm::vec3(0.0f, 1.0f, 0.0f);
GLfloat yawStart = 0.0f;
GLfloat pitchStart = 0.0f;
GLdouble lastXStart = (GLdouble)WIDTH / 2.0;
GLdouble lastYStart = (GLdouble)HEIGHT / 2.0;
GLfloat fovStart = 50.0f;

// Camera vars
glm::vec3 cameraPos = cameraPosStart;
glm::vec3 cameraFront = cameraFrontStart;
glm::vec3 cameraUp = cameraUpStart;
GLfloat yaw = yawStart;
GLfloat pitch = pitchStart;
GLdouble lastX = lastXStart;
GLdouble lastY = lastYStart;
GLfloat fov = fovStart;
bool keys[1024];
GLfloat currIsoVal = 40.0f;

// Deltatime - normalizes between systems
GLfloat currentFrame = 0.0f;
GLfloat deltaTime = 0.0f;	// Time between current frame and last frame
GLfloat lastFrame = 0.0f;  	// Time of last frame
GLfloat saveLast = 0.0f;

// Function prototypes
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, GLdouble xoffset, GLdouble yoffset);
void mouse_callback(GLFWwindow* window, GLdouble xpos, GLdouble ypos);
void handleGLerrors(GLenum error);
void do_movement();
void reset();
std::string findBrickBinaryFile(int brick_number, int timestep);
void fillVolumeBrickMap();

// Shader files
const char* vertFile = ".\\raycaster.vert";
const char* fragFile = ".\\raycaster.frag";
//const char* vertFile = ".\\isosurface.vert";
//const char* fragFile = ".\\isosurface.frag";
//const char* fragFile = ".\\isosurface.frag";
ShaderProgram shaderProgram;

int64_t BRICK_DIM[3] = { 256, 256, 128 };
int64_t VOLUME_GRID[3] = { 8, 8, 15 };
  
std::map<int, dim3D> volumeBrickMap;


int main() {
	
	/****** Set up windowing & UI with GLFW ******/

	// Instantiate the GLFW window
	glfwInit();

	// Set all required options for GLFW
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	glfwWindowHint(GLFW_DOUBLEBUFFER, GL_TRUE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // for mac

	// Create a window object
	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Volume Visualization", nullptr, nullptr);
	if (window == nullptr)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Set the required callback functions
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

	// Options
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);



	/***** Initialize GLEW *******/

	// Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
	glewExperimental = GL_TRUE;

	// Initialize GLEW to setup the OpenGL Function pointers
	if (glewInit() != GLEW_OK)
	{
		std::cout << "Failed to initialize GLEW" << std::endl;
		std::cin.ignore();
		return -1;
	}

	// Define viewpoint dimensions
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height); // location of lower left corner of window



	/***** Store Volume Data as 3D Texture *****/

	fillVolumeBrickMap();
	VolumeReader<GLubyte> volume(BRICK_DIM, VOLUME_GRID,
		findBrickBinaryFile, &volumeBrickMap);
	
	// Load datasets to GPU - No Compression
	//bool texLoadSuccess = volume.LoadBricksToTexture(448, 8, 8, 7, 273, true, true); // max memory

	// Compress
	//bool texLoadSuccess = volume.LoadBricksToTexture(256, 8, 8, 4, 273, false);
	//bool texLoadSuccess = volume.LoadBrickToTexture(700, 273, false, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(384, 8, 8, 6, 273, false, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(256, 8, 8, 4, 273, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(448, 8, 8, 7, 270, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(960, 8, 8, 15, 270, false);

	//VolumeKdtree * myTree = new VolumeKdtree(volume.data, volume.dataDims[0], volume.dataDims[1], volume.dataDims[2]);

	//MidRangeTree * myTree = new MidRangeTree();
	//myTree->open("mid_range_tree_1brick.bin");
	//myTree->open("mid_range_tree_256bricks.bin")
	//myTree->open("tree_384_6tolerance.bin");

	
	//myTree->setMaxEpochs(1);
	//myTree->setErrorTolerance(4);
	//myTree->setErrorTolerance(6);

	//DebugTimer::Begin(1, "TOTAL_CONSTRUCTION");
	//myTree->build(true);
	//DebugTimer::End("TOTAL_CONSTRUCTION");
	
	

	//myTree->save("tree_256brick_2error.bin");
	//myTree->save("mid_range_tree_256bricks.bin");
	//myTree->save("mid_range_tree_256bricks_2error.bin");
	//myTree->save("mid_range_tree_448bricks_0error.bin");

	/*
	std::vector<unsigned char> treeData;
	myTree->levelCut(myTree->maxTreeDepth, treeData);
	//delete myTree; 

	
	// Create OpenGL 3D Texture
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_3D, textureId);

	// Set the texture parameters
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Load uncompressed data vector into 3D Texture
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, (GLsizei)BRICK_DIM[0], (GLsizei)BRICK_DIM[1], (GLsizei)BRICK_DIM[2], 0, GL_RED, GL_UNSIGNED_BYTE, &treeData[0]);
	//glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, (GLsizei)volume.dataDims[0], (GLsizei)volume.dataDims[1], (GLsizei)volume.dataDims[2], 0, GL_RED, GL_UNSIGNED_BYTE, &treeData[0]);
	*/


	/*
	std::vector<unsigned char> compressedVolume;
	myTree->convertToByteArray(compressedVolume);   
	std::cout << compressedVolume.size() << std::endl;

	// Store the compressed tree into 1D Texture
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_1D, textureId);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_R, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_R, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	
	
	glTexImage1D(GL_TEXTURE_1D, 0, GL_INTENSITY, compressedVolume.size(), 0, GL_RED, GL_UNSIGNED_BYTE, &compressedVolume[0]);
	*/

	// Store the distance maps as uniform buffer objects
	//const int MAX_HEIGHT = 40;
	//GLuint myArrayUBO;
	//glGenBuffers(1, &myArrayUBO);
	//glBindBuffer(GL_UNIFORM_BUFFER, myArrayUBO);
	//glBufferData(GL_UNIFORM_BUFFER, sizeof(GLint) * MAX_HEIGHT, &myTree->distanceMap[0], GL_DYNAMIC_DRAW);

	//GLuint myArrayUBO1;
	//glGenBuffers(1, &myArrayUBO1);
	//glBindBuffer(GL_UNIFORM_BUFFER, myArrayUBO1);
	//glBufferData(GL_UNIFORM_BUFFER, sizeof(GLint) * MAX_HEIGHT, &myTree->distanceMap_range[0], GL_DYNAMIC_DRAW);

	//bool texLoadSuccess = volume.LoadBricksToTexture(512, 8, 8, 8, 273, true);
	//bool texLoadSuccess = volume.LoadBrickToTexture(700, 273, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(256, 8, 8, 4, 273, false);
	bool texLoadSuccess = volume.LoadBricksToTexture(384, 8, 8, 6, 273, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(448, 8, 8, 7, 273, false);
	//bool texLoadSuccess = volume.LoadBricksToTexture(960, 8, 8, 15, 273, false, false);
	//std::cout << volume.dataDims[0] << " " << volume.dataDims[1] << " " << volume.dataDims[2] << std::endl;
	//std::cout << volume.data.size() << std::endl;
	// test //
	//VolumeKdtree * myTree = new VolumeKdtree();

	
	VolumeKdtree * myTree = new VolumeKdtree(volume.data, volume.dataDims[0], volume.dataDims[1], volume.dataDims[2]);
	//MidRangeTree * myTree = new MidRangeTree(volume.data, volume.dataDims[0], volume.dataDims[1], volume.dataDims[2]);
	myTree->setMaxEpochs(2);
	myTree->setErrorTolerance(1);
	
	
    DebugTimer::Begin(1, "TOTAL_CONSTRUCTION");
	myTree->build(true);
	DebugTimer::End("TOTAL_CONSTRUCTION");
	
	

	//myTree->save("mid_range_tree_allbricks.bin");

	//myTree->open("tree_384_2tolerance_2epoch.bin");
	//myTree->open("tree_384_6tolerance_recovered.bin");
	myTree->save("tree_384_1tolerance.bin");
	//myTree->save("tree_256_12tolerance_r.bin");
	//myTree->open("tree_brick_5extralevels.bin");
	//myTree->save("tree_brick_preorder_2extralevels.bin");
	//myTree->save("tree_384_preorder_0extralevels.bin");
	//myTree->save("tree256_hashed.bin");
	//myTree->save("tree_brick_hashed.bin");
	//myTree->open("tree_brick_hashed.bin");
	//myTree->open("tree_brick.bin");
	//myTree->save("tree256.bin");
	//myTree->open("tree256.bin");
	//myTree->open("tree_brick_preorder.bin");

	std::vector<unsigned char> treeData;
	myTree->levelCut(myTree->maxTreeDepth, treeData);

	//std::cout << "MAX ERROR: " << myTree->measureMaxError() << std::endl;
	//std::cout << "MEAN ERROR: " << myTree->measureMeanError() << std::endl;
	//std::vector<unsigned char> errorData;
	//myTree->queryError(errorData);
	//std::cout << treeData.size() << " " << volume.data.size() << std::endl;
	//glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, (GLsizei)volume.dataDims[0], (GLsizei)volume.dataDims[1], (GLsizei)volume.dataDims[2], 0, GL_RED, GL_UNSIGNED_BYTE, &volume.data[0]);

	glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, (GLsizei)volume.dataDims[0], (GLsizei)volume.dataDims[1], (GLsizei)volume.dataDims[2], 0, GL_RED, GL_UNSIGNED_BYTE, &treeData[0]);
	
	//glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, (GLsizei)volume.dataDims[0], (GLsizei)volume.dataDims[1], (GLsizei)volume.dataDims[2], 0, GL_RED, GL_UNSIGNED_BYTE, &errorData[0]);
	//std::vector<unsigned char> diff;
	//std::set_difference(volume.data.begin(), volume.data.end(), treeData.begin(), treeData.end(),
	//	std::inserter(diff, diff.begin()));
	//int64_t sumDiff = std::accumulate(diff.begin(), diff.end(), 0);
	//std::cout << " ERROR = " << sumDiff << std::endl;

	//delete myTree;
	

	//if (!texLoadSuccess) {
	//	std::cin.ignore();
	//	exit(-1);
	//}

	//std::cout << "texture check" << std::endl;
	handleGLerrors(glGetError());


	/****** Set up Shaders ******/

	
	shaderProgram.loadShaderFromFile(vertFile, GL_VERTEX_SHADER);
	shaderProgram.loadShaderFromFile(fragFile, GL_FRAGMENT_SHADER);
	shaderProgram.createProgram();
	shaderProgram.Use();

	shaderProgram.addUniform("ViewMatrix");
	shaderProgram.addUniform("ProjectionMatrix");
	shaderProgram.addUniform("volume");
	shaderProgram.addUniform("camPos");
	shaderProgram.addUniform("step_size");
	shaderProgram.addUniform("TransformationMatrix");
	shaderProgram.addUniform("isoValue");
	shaderProgram.addUniform("origTreeDepth");
	shaderProgram.addUniform("maxTreeDepth");

	// Note: step size should be adjusted based on volume size, but be careful of fps
	glUniform3f(shaderProgram.uniforms["step_size"], 1.0f / (float)BRICK_DIM[0],
		1.0f / (float)BRICK_DIM[1], 1.0f / (float)BRICK_DIM[2]);
	glUniform1i(shaderProgram.uniforms["volume"], 0);
	
	glUniform1f(shaderProgram.uniforms["isoValue"], currIsoVal / 255.0f);
	//glUniform1i(shaderProgram.uniforms["origTreeDepth"], myTree->origTreeDepth);
	//glUniform1i(shaderProgram.uniforms["maxTreeDepth"], myTree->maxTreeDepth);

	glm::mat4 trans;
	//trans = glm::rotate(trans, glm::radians(50.0f), glm::vec3(1.0, 0.0, 0.0));
	//trans = glm::scale(trans, glm::vec3(2.0, 2.0, 2.0));
	glUniformMatrix4fv(shaderProgram.uniforms["TransformationMatrix"], 1, GL_FALSE, glm::value_ptr(trans));

	//GLuint distanceMapIdx = glGetUniformBlockIndex(shaderProgram.programId, "distanceMap");
	//glUniformBlockBinding(shaderProgram.programId, distanceMapIdx, 0);
	//glBindBufferBase(GL_UNIFORM_BUFFER, 0, myArrayUBO);

	//GLuint distanceMapRangeIdx = glGetUniformBlockIndex(shaderProgram.programId, "distanceMapRange");
	//glUniformBlockBinding(shaderProgram.programId, distanceMapRangeIdx, 1);
	//glBindBufferBase(GL_UNIFORM_BUFFER, 1, myArrayUBO1);

	//std::cout << "shader check" << std::endl;
	handleGLerrors(glGetError());



	/**** Set up VAO & VBO to render a unit cube *****/

	UnitBrick brick;
	brick.Setup();

	//std::cout << "vertex check" << std::endl;
	handleGLerrors(glGetError());
	
	
	/***** Render in Game Loop ******/

	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	brick.Bind();

	while (!glfwWindowShouldClose(window))
	{
		// Time rendering
		glfwSwapInterval(0); //vertical sync
		glFinish();
		DebugTimer::Begin(500, "LOOP");

		// Calculate deltatime of current frame
		currentFrame = (float)glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		saveLast = lastFrame;
		lastFrame = currentFrame;

		// check if any events have been activated
		glfwPollEvents();
		do_movement();

		// Clear color & depth buffer
		//glClearColor(0.0f, 0.0f, 0.0f, 1.0f); //black
		glClearColor(255.0f, 255.0f, 255.0f, 1.0f); //white
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// dynamic matrices
		glm::mat4 viewMatrix = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		glm::mat4 projMatrix = glm::perspectiveFov(glm::radians(fov), (GLfloat)width, (GLfloat)height, 0.1f, 100.f);

		// Update uniform data
		glUniformMatrix4fv(shaderProgram.uniforms["ProjectionMatrix"], 1, GL_FALSE, glm::value_ptr(projMatrix));
		glUniformMatrix4fv(shaderProgram.uniforms["ViewMatrix"], 1, GL_FALSE, glm::value_ptr(viewMatrix));
		glUniform3fv(shaderProgram.uniforms["camPos"], 1, &(cameraPos.x));

		brick.Draw();

		handleGLerrors(glGetError());

		glFinish();
    		DebugTimer::End("LOOP");
		glfwSwapBuffers(window);
	}

	// De-allocate resources
	shaderProgram.StopUse();
	brick.Unbind();
	brick.Delete();

	// Clean resources
	glfwTerminate();
	
	return 0;
}

void handleGLerrors(GLenum error) {

	if (error == GL_NO_ERROR) {
	}
	else if (error == GL_INVALID_ENUM) {
		std::cout << "invalid enum" << std::endl;
	}
	else if (error == GL_INVALID_VALUE) {
		std::cout << "invalid value" << std::endl;
	}
	else if (error == GL_INVALID_OPERATION) {
		std::cout << "invalid operation" << std::endl;
	}
	else if (error == GL_STACK_OVERFLOW) {
		std::cout << "STACK OVERFLOW" << std::endl;
	}
	else if (error == GL_STACK_UNDERFLOW) {
		std::cout << "STACK UNDERFLOW" << std::endl;
	}
	else if (error == GL_OUT_OF_MEMORY) {
		std::cout << "out of memory" << std::endl;
	}
	else if (error == GL_INVALID_FRAMEBUFFER_OPERATION) {
		std::cout << "invalid frame buffer operation" << std::endl;
	}
	else if (error == GL_CONTEXT_LOST) {
		std::cout << "context lost" << std::endl;
	}
	else if (error == GL_TABLE_TOO_LARGE) {
		std::cout << "table too large" << std::endl;
	}
	else {
		std::cout << "unhandled GL error" << std::endl;
	}

}

// Moves/alters the camera positions based on user input
void do_movement()
{
	// Camera controls
	//GLfloat cameraSpeed = 0.08f;
	GLfloat cameraSpeed = 2.5f * deltaTime;
	if (keys[GLFW_KEY_UP])
		cameraPos += cameraSpeed * cameraFront;
	if (keys[GLFW_KEY_DOWN])
		cameraPos -= cameraSpeed * cameraFront;
	if (keys[GLFW_KEY_LEFT])
		cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	if (keys[GLFW_KEY_RIGHT])
		cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	
	
		
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	if (key == GLFW_KEY_ENTER && action == GLFW_PRESS)
		reset();

	if (key == GLFW_KEY_0 && action == GLFW_PRESS) {
		currIsoVal = std::max(0.0f, currIsoVal - 5.0f);
		glUniform1f(shaderProgram.uniforms["isoValue"], currIsoVal / 255.0f);
		std::cout << "Isovalue: " << currIsoVal << std::endl;
	}
	if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
		currIsoVal = std::min(255.0f, currIsoVal + 5.0f);
		glUniform1f(shaderProgram.uniforms["isoValue"], currIsoVal / 255.0f);
		std::cout << "Isovalue: " << currIsoVal << std::endl;
	}

	if (key >= 0 && key < 1024)
	{
		if (action == GLFW_PRESS)
			keys[key] = true;
		else if (action == GLFW_RELEASE)
			keys[key] = false;
	}
}

void scroll_callback(GLFWwindow* window, GLdouble xoffset, GLdouble yoffset)
{

	if (fov >= 1.0f && fov <= fovStart)
		fov -= (GLfloat)yoffset;
	if (fov <= 1.0f)
		fov = 1.0f;
	if (fov >= fovStart)
		fov = fovStart;
}


/*
Look around in direction of mouse when mouse button is clicked.
*/
bool firstMouse = true;
void mouse_callback(GLFWwindow* window, GLdouble xpos, GLdouble ypos)
{
	GLint mb = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1);

	if (mb) {
		if (firstMouse)
		{
			lastX = xpos;
			lastY = ypos;
			firstMouse = false;
		}

		GLdouble xoffset = xpos - lastX;
		GLdouble yoffset = lastY - ypos; // Reversed since y-coordinates go from bottom to left
		lastX = xpos;
		lastY = ypos;

		GLfloat sensitivity = 1.0;	// adjust
		xoffset *= sensitivity;
		yoffset *= sensitivity;

		pitch += (GLfloat)yoffset;
		yaw += (GLfloat)xoffset; //neg

								 // Make sure when pitch is out of bounds, screen doesn't get flipped
		if (pitch > 89.0f)
			pitch = 89.0f;
		if (pitch < -89.0f)
			pitch = -89.0f;

		glm::vec3 front;
		front.x = cos(glm::radians(pitch)) * cos(glm::radians(yaw));
		front.y = sin(glm::radians(pitch));
		front.z = sin(glm::radians(yaw));
		cameraFront = glm::normalize(front);

	}
	else {
		firstMouse = true;
	}

}

void reset()
{
	cameraPos = cameraPosStart;
	cameraFront = cameraFrontStart;
	cameraUp = cameraUpStart;
	yaw = yawStart;
	pitch = pitchStart;
	lastX = lastXStart;
	lastY = lastYStart;
	fov = fovStart;
}

std::string findBrickBinaryFile(int brick_number, int timestep) {

	// string manipulations 
	std::stringstream output;
	std::stringstream timestep3;
	std::stringstream timestep4;
	std::stringstream brick4;

	timestep3 << std::setw(3) << std::setfill('0') << timestep;
	timestep4 << std::setw(4) << std::setfill('0') << timestep;
	brick4 << std::setw(4) << std::setfill('0') << brick_number;

	std::string topDir = "C:\\Users\\tlarrue\\Documents\\richtmyer-meshkov_instability_dataset\\all_bricks\\";

	output << topDir << "bob" << timestep3.str() << "\\d_" << timestep4.str() << "_" << brick4.str();

	return output.str();
}

void fillVolumeBrickMap() {
	int i, j, k;
	i = j = k = 0;
	for (int b = 0; b < 8 * 8 * 15; b++) {
		volumeBrickMap[b][0] = i;
		volumeBrickMap[b][1] = j;
		volumeBrickMap[b][2] = k;
		
		if ((b + 1) % 64 == 0) {
			i = j = 0;
			k++;
		}
		else if ((b + 1) % 8 == 0) {
			i = 0;
			j++;
		}
		else {
			i++;
		}
	}
}