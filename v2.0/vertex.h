#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include "edge.h"
#include "primitive.h"
#include "plane.h"
struct Edge;
struct Primitive;

// Per-vertex data
struct Vertex {
	int id = -1;
	glm::dvec3 pos;
	std::vector<Primitive*> triangles;
	std::vector<Edge*> edges;
	bool silhouette = false;
	bool splitline = false;


	bool isInFrontOfVertex(Vertex* v, Plane* plane) const {
		if (v->id == id) return false;
		return fabs(plane->distToPoint(pos)) < fabs(plane->distToPoint(v->pos));
	}

	struct cmp_ptr
	{
		bool operator()(const Vertex* lhs, const Vertex* rhs) const
		{
			return lhs->id < rhs->id;
		}
	};
};

struct VertexVis {
	glm::dvec3 pos;
	glm::dvec3 normal;
	glm::vec3 color = glm::vec3(1);
	float tri_id;
	//float visible = 1.f;
};