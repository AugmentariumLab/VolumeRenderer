#include <vector>
#include <map>

#include "VolumeReader.h"
#include "ShaderProgram.h"
#include "UnitBrick.h"
#include "DebugTimer.h"

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
glm::vec3 cameraPosStart = glm::vec3(0.0f, 0.0f, -5.0f);
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

unsigned int VOLUME_DIM[3] = { 256, 256, 128 };
unsigned int VOLUME_GRID[3] = { 8, 8, 15 };

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
	VolumeReader<GLubyte> volume(VOLUME_DIM, VOLUME_GRID,
		findBrickBinaryFile, &volumeBrickMap);

	//bool texLoadSuccess = volume.LoadBrickToTexture(700, 273);
	bool texLoadSuccess = volume.LoadBricksToTexture(448, 8, 8, 7, 273);

	if (!texLoadSuccess) {
		std::cin.ignore();
		exit(-1);
	}

	//std::cout << "texture check" << std::endl;
	handleGLerrors(glGetError());


	/****** Set up Shaders ******/

	ShaderProgram shaderProgram;
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

	glUniform3f(shaderProgram.uniforms["step_size"], 1.0f / (float)VOLUME_DIM[0], 
		1.0f / (float)VOLUME_DIM[1], 1.0f / (float)VOLUME_DIM[2]);
	glUniform1i(shaderProgram.uniforms["volume"], 0);


	glm::mat4 trans;
	trans = glm::scale(trans, glm::vec3(2.0, 2.0, 2.0));
	glUniformMatrix4fv(shaderProgram.uniforms["TransformationMatrix"], 1, GL_FALSE, glm::value_ptr(trans));

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
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f); //black
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
	GLfloat cameraSpeed = 5.0f * deltaTime;
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