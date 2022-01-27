#version 410 core

// Output for on-screen color
layout(location = 0) out float outColor;

// Interpolated output data from vertex shader
in float fragId;

void main() 
{
	outColor = fragId;
}
