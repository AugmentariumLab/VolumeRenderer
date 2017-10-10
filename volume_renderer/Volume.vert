#version 330 core


uniform mat4 ProjectionMatrix, ViewMatrix;
//uniform mat4 ProjectionMatrix, ViewMatrix, ModelMatrix, transform;
uniform mat4 MVP;
layout (location = 0) in vec3 position;

void main()
{
	gl_Position = ProjectionMatrix * ViewMatrix * vec4(position, 1.f);
	//gl_Position = ProjectionMatrix * ViewMatrix * ModelMatrix * transform * vec4(position, 1.f);
	//gl_Position = MVP * vec4(position, 1.f);
}