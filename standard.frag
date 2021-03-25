#version 330 core

// uniform mat4 lightMVP;
uniform vec3 lightPos = vec3(-3,3,-3);
uniform vec3 lightPos2 = vec3(-3,3,3);
uniform vec3 lightPos3 = vec3(3,-3,0);
// Output for on-screen color
layout(location = 0) out vec4 outColor;

// Interpolated output data from vertex shader
in vec3 fragPos;    // World-space position
in vec3 fragNormal; // World-space normal
in float fragSelect;
//in float id;

void main() {

	vec3 lightDir = normalize(lightPos - fragPos);
	vec3 lightDir2 = normalize(lightPos2 - fragPos);
	vec3 lightDir3 = normalize(lightPos3 - fragPos);

	vec3 color = vec3(1.f, 1.f, 1.f);
	if (fragSelect > 1E-6) color = vec3(1.f, 0.f, 0.f);
	vec3 diffuse = vec3(max(dot(fragNormal, lightDir), 0.0));
	vec3 diffuse2 = vec3(max(dot(fragNormal, lightDir2), 0.0));
	vec3 diffuse3 = vec3(max(dot(fragNormal, lightDir3), 0.0));
	outColor = vec4((color*(0.8f*(diffuse+diffuse2+diffuse3))), 1.0);
}