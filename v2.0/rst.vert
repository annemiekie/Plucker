#version 410 core

// Model/view/projection matrix
uniform mat4 mvp;
uniform int secondpass;

// Per-vertex attributes
layout(location = 0) in vec3 pos; // World-space position
layout(location = 1) in vec3 normal; // World-space normal
layout(location = 2) in vec3 bary;
layout(location = 3) in vec3 center;
layout(location = 5) in float id; // prim id

// Data to pass to fragment shader
out float fragId;
out vec3 fragBary;

void main() {
	// Transform 3D position into on-screen position
    gl_Position = mvp * vec4(pos, 1.0);

    // Pass position and normal through to fragment shader
    fragBary = bary;
    fragId = id;
}
   // if (secondpass==1) newpos = center + 1.1f*(pos-center);
