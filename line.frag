#version 330 core

// Output for on-screen color
layout(location = 0) out vec4 outColor;

// Interpolated output data from vertex shader
in vec4 lineColor; // Color

void main() 
{
	outColor = lineColor;
}
