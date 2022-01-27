#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>
struct Edge;

// Per-vertex data
struct Vertex {
	glm::vec3 pos;
	int id=-1;
	std::vector<int> triangles;
	std::vector<Edge*> edges;

	struct cmp_ptr
	{
		bool operator()(const Vertex* lhs, const Vertex* rhs) const
		{
			return lhs->id < rhs->id;
		}
	};
};

struct VertexVis {
	glm::vec3 pos;
	glm::vec3 normal;
	glm::vec3 color;
	//glm::vec3 center;
	float tri_id;
};