#version 330 core

// Output for on-screen color
layout(location = 0) out vec4 outColor;

// Interpolated output data from vertex shader
in vec3 lineColor; // Color

void main() 
{
	outColor = vec4(lineColor,1);
}
