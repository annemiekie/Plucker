#version 330 core

// Model/view/projection matrix
uniform mat4 mvp;
uniform float alpha;
uniform int setcol;
uniform vec3 setcolor;

// Per-vertex attributes
layout(location = 0) in vec3 pos; // World-space position
layout(location = 1) in vec3 color; // Color

// Data to pass to fragment shader
out vec4 lineColor;

void main() {
	// Transform 3D position into on-screen position
    gl_Position = mvp * vec4(pos, 1.0);
    vec4 lc;
    if (setcol==1) lc = vec4(setcolor, 1.0);
    else lc = vec4(color, alpha);
    lineColor = lc;
}
