#version 410 core

// Output for on-screen color
layout(location = 0) out float outColor;

// Interpolated output data from vertex shader
in float fragId;

void main() 
{
	int r = ((int)fragId & 0x000000FF) >>  0;
	int g = ((int)fragId & 0x0000FF00) >>  8;
	int b = ((int)fragId & 0x00FF0000) >> 16;
	outColor = vec4(r/255.0f, g/255.0f, b/255.0f,1.0);
}