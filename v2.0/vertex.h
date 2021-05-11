#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>

// Per-vertex data
struct Vertex {
	glm::vec3 pos;
	glm::vec3 normal;
	glm::vec3 bary;
	glm::vec3 center;
	float selected;
	float id;
};