#version 330 core

// uniform mat4 lightMVP;
uniform vec3 lightPos = vec3(-3,3,-3);

// Output for on-screen color
layout(location = 0) out vec4 outColor;

// Interpolated output data from vertex shader
in vec3 fragPos;    // World-space position
in vec3 fragNormal; // World-space normal

void main() {

	vec3 lightDir = normalize(lightPos - fragPos);
	outColor = vec4(vec3(max(dot(fragNormal, lightDir), 0.0)), 1.0); 
}