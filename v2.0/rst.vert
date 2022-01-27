#version 410 core

// Model/view/projection matrix
uniform mat4 mvp;

// Per-vertex attributes
layout(location = 0) in vec3 pos; // World-space position
layout(location = 3) in float id; // prim id

// Data to pass to fragment shader
out float fragId;

void main() {
	// Transform 3D position into on-screen position
    gl_Position = mvp * vec4(pos, 1.0);
    fragId = id;
}
