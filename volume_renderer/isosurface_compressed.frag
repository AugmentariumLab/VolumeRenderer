#version 430 core

layout(location = 0) out vec4 vFragColor;	//fragment shader output

smooth in vec3 vUV;				//3D texture coordinates form vertex shader 
								//interpolated by rasterizer

//uniforms
uniform vec3		camPos;			//camera position
uniform float		isoValue;		// dynamic isovalue

//constants
const int MAX_SAMPLES = 300;		//total samples for each ray march step
const vec3 texMin = vec3(0);		//minimum texture access coordinate
const vec3 texMax = vec3(1);		//maximum texture access coordinate
const float DELTA = 0.01;			//the step size for gradient calculation

 layout(std430, binding = 3) buffer layoutName
 {
     int data_SSBO[];
 };

struct node {
	double idx;
	vec3 boundMin;
	vec3 boundMax;
	uint scalarMin;
	uint scalarMax;
};

void main (void) {
  //get the 3D texture coordinates for lookup into the volume dataset
  vec3 dataPos = vUV;

  //Gettting the ray marching direction:
  //get the object space position by subracting 0.5 from the
  //3D texture coordinates. Then subtraact it from camera position
  //and normalize to get the ray marching direction
  vec3 geomDir = normalize((vUV-vec3(0.5)) - camPos); 

  // Start at root
  //node currentNode = node(0, vec3(0,0,0), vec3(1,1,1), )

  vFragColor = vec4(0.5, 0.5, 0.5, 0.5);
}