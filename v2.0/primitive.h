#ifndef PRIMITIVE_H
#define PRIMITIVE_H
#include <vector>
#include "vertex.h"
#include "ray.h"
#include "plane.h"

struct Edge;

struct Primitive {
	int id;
	glm::dvec3 normal;
	glm::dvec3 center;
	Vertex* vertices[3];
	Edge* edges[3];
	Ray rays[3];

	float getIntersectionDepth(const Ray& r) {
		double ndotdir = glm::dot(normal, r.direction);
		return glm::dot(normal, vertices[0]->pos - r.origin) / ndotdir;
	}

	bool onRightSideOfSplitEdge(Ray& split, bool side) {
		Ray checkRay = Ray(center + normal, center);
		return split.side(checkRay) == side;
	}

	bool onRightSideOfSplitVertex(Vertex* v, Ray& split , bool side) {
		// check if center of primitive lies on correct side of splitline edge
		if (split.side(Ray(center + normal, center)) == side) return true;

		// check if one of vertices lies on correct side of splitline edge
		for (Vertex* v2 : vertices)
			if (v != v2 && split.side(Ray(v2->pos + normal, v2->pos)) == side)
				return true;
		return false;
	}

	bool intersection(const Ray& or, float& depth) {
		for (Ray& r : rays) if (r.side(or)) return false;
		depth = getIntersectionDepth(or);
		return true;
	}

	Plane getPlane() {
		return Plane(vertices[0]->pos, normal);
	}

	std::vector<Ray> getRayVector() {
		return  { rays[0], rays[1], rays[2] };
	}

	bool hasVertex(int id) {
		for (Vertex* pv : vertices) {
			if (pv->id == id) return true;
		}
		return false;
	}

	bool vertexOnPositiveSidePlane(Vertex* v) {
		Plane plane = getPlane();
		return plane.pointOnPositiveSide(v->pos) && !plane.pointOnPlane(v->pos, 1E-7) && !hasVertex(v->id);
	}
	//bool vertexInCommonWithEdge(Edge* e) {
	//	for (Vertex* pv : vertices) {
	//		for (Vertex* ev : e->vertices) {
	//			if (pv->id == ev->id) return true;
	//		}
	//	}
	//	return false;
	//}
};

#endif
