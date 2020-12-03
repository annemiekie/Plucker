#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>

// Per-vertex data
struct Vertex {
	glm::vec3 pos;
	glm::vec3 normal;
	float id;
};