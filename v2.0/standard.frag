#version 330 core

// uniform mat4 lightMVP;
uniform vec3 lightPos1;
uniform vec3 lightPos2;
uniform vec3 lightPos3;
// Output for on-screen color
layout(location = 0) out vec4 outColor;

// Interpolated output data from vertex shader
in vec3 fragPos; 
in vec3 fragNormal;
in vec3 fragColor;

void main() {

	vec3 lightDir1 = normalize(lightPos1 - fragPos);
	vec3 lightDir2 = normalize(lightPos2 - fragPos);
	vec3 lightDir3 = normalize(lightPos3 - fragPos);

	vec3 diffuse1 = vec3(max(dot(fragNormal, lightDir1), 0.0));
	vec3 diffuse2 = vec3(max(dot(fragNormal, lightDir2), 0.0));
	vec3 diffuse3 = vec3(max(dot(fragNormal, lightDir3), 0.0));
	outColor = vec4((fragColor*(0.5f*(diffuse1+diffuse2+diffuse3))),1.0);
}