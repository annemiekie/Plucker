#ifndef PRIMITIVE_H
#define PRIMITIVE_H
#include <vector>
#include <set>
#include "edge.h"
#include "vertex.h"

struct Primitive {
	int index;
	glm::vec3 normal;
	std::vector<Vertex> vertices;
	mutable std::vector<Edge> edges = std::vector<Edge>();
};

#endif
