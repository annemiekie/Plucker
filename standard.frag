#version 330 core

// Global variables for lighting calculations
uniform vec3 viewPos;

// uniform mat4 lightMVP;
uniform vec3 lightPos = vec3(-3,-3,-3);

// Output for on-screen color
layout(location = 0) out vec4 outColor;

// Interpolated output data from vertex shader
in vec3 fragPos;    // World-space position
in vec3 fragNormal; // World-space normal
in int gl_PrimitiveID;
in vec4 includes;

void main() {

	//if (gl_PrimitiveID == includes.x || gl_PrimitiveID == includes.y || gl_PrimitiveID == includes.z || gl_PrimitiveID == includes.w) {
	// Output the normal as color
		vec3 lightDir = normalize(lightPos - fragPos);
		outColor = vec4(vec3(max(dot(fragNormal, lightDir), 0.0)), 1.0); 
	//}
	//outColor = vec4(1.0, 0.0, 0.0, 1.0); 

}