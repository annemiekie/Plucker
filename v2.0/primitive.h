#ifndef PRIMITIVE_H
#define PRIMITIVE_H
#include <vector>
//#include "edge.h"
#include "vertex.h"
#include "ray.h"
#include "plane.h"

struct Edge;

struct Primitive {
	int id;
	glm::vec3 normal;
	glm::vec3 center;
	Vertex* vertices[3];
	Edge* edges[3];
	Ray rays[3];

	float getIntersectionDepth(const Ray& r) {
		double ndotdir = glm::dot(normal, glm::vec3(r.direction));
		return glm::dot(normal, vertices[0]->pos - glm::vec3(r.origin)) / ndotdir;
	}

	Plane getPlane() {
		return Plane(vertices[0]->pos, normal);
	}
};

#endif
