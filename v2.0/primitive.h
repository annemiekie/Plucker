#ifndef PRIMITIVE_H
#define PRIMITIVE_H
#include <vector>
#include "vertex.h"
#include "ray.h"
#include "plane.h"
#include <set>
struct Edge;

struct Primitive {
	int id;
	glm::dvec3 normal;
	glm::dvec3 center;
	Vertex* vertices[3];
	Edge* edges[3];
	Ray rays[3];
	Plane* plane = NULL;
	std::set<int> esls;

	double getIntersectionDepth(const Ray& r, bool inplane = false, bool closestEdge = false) {
		//if (inplane || getPlane().rayInPlane(r, 1E-7)) {
		//	double depth = -INFINITY;
		//	for (Ray er : rays) {
		//		double edgeIntersect = r.depthToIntersectionWithRay(er);
		//		if (closestEdge) {
		//			if (edgeIntersect < depth) depth = edgeIntersect;
		//		}
		//		else if (edgeIntersect > depth) depth = edgeIntersect;
		//	}
		//	return depth;
		//}
		double ndotdir = glm::dot(normal, r.direction);
		return glm::dot(normal, vertices[0]->pos - r.origin) / ndotdir;
	}

	//double getIntersectionDepth(const Ray& r, bool inplane = false, bool closestEdge = false) {
	//	return getPlane().rayIntersectionDepth(r);
	//}
	bool vertexNeighbor(Primitive* other) {
		for (Vertex* v1 : vertices) {
			for (Vertex* v2 : other->vertices) {
				if (v1->id == v2->id) return true;
			}
		}
		return false;
	}

	std::set<Vertex*> sameVertices(Primitive* other) {
		std::set<Vertex*> verts;
		for (Vertex* v1 : vertices) {
			for (Vertex* v2 : other->vertices) {
				if (v1->id == v2->id) verts.insert(v1);
			}
		}
		return verts;
	}

	//bool edgeNeighbor(Primitive* other) {
	//	for (auto* e1 : edges) {
	//		for (auto* e2 : other->edges) {
	//			if (e1->id == e2->id) return true;
	//		}
	//	}
	//	return false;
	//}

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

	bool intersection(const Ray& or , double& depth, bool inplane = false) {
		for (Ray& r : rays) if (r.side(or)) return false;
		depth = getIntersectionDepth(or, inplane);
		return true;
	}

	Plane* getPlane() {
		if (plane == NULL) plane = new Plane(vertices[0]->pos, normal);
		return plane;
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
