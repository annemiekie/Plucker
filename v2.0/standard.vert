#version 330 core

// Model/view/projection matrix
uniform mat4 mvp;

// Per-vertex attributes
layout(location = 0) in vec3 pos; 
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 color;

// Data to pass to fragment shader
out vec3 fragPos;
out vec3 fragNormal;
out vec3 fragColor;

void main() {
	// Transform 3D position into on-screen position
    gl_Position = mvp * vec4(pos, 1.0);

    // Pass position and normal through to fragment shader
    fragPos = pos;
    fragNormal = normal;
    fragColor = color;

    gl_PointSize = 2;
}
