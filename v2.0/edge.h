#ifndef EDGE_H
#define EDGE_H
#include "plane.h"
#include "primitive.h"
#include "vertex.h"

struct Edge {
	mutable int id = -1;
	Vertex* vertices[2];
	Ray ray;
	mutable std::vector<Primitive*> triangles;
	bool silhouette = false;
	bool splitline = false;

	bool intersectsRay(Ray& oray) {
		glm::dvec3 ipoint = ray.pointOfintersectWithRay(oray);
		glm::dvec3 tvec = (ipoint - vertices[0]->pos) /
						 (vertices[1]->pos - vertices[0]->pos);
		return (tvec.x > 0 && tvec.x < 1);
	}

	bool pointOnRay(glm::dvec3& pt) {
		glm::vec3 tvec = (pt - vertices[0]->pos) /
						(vertices[1]->pos - vertices[0]->pos);
		return (tvec.x > 0 && tvec.x < 1);
	}

	bool isSilhouetteForPos(glm::dvec3& pt) const {
		bool side;
		return isSilhouetteForPos(pt, side);
	}

	bool isSilhouetteForPos(glm::dvec3& pt, bool& side) const {
		if (triangles.size() == 1) return true;
		Plane p1 = triangles[0]->getPlane();
		Plane p2 = triangles[1]->getPlane();
		return isSilhouetteForPos(pt, p1, p2, side);
	}

	bool isSilhouetteForPos(glm::dvec3& pt, Plane& p1, Plane& p2, bool& side) const {
		side = p1.pointOnPositiveSide(pt);// || p1.pointOnPlane(pt, 1E-6);
		bool side2 = p2.pointOnPositiveSide(pt);// || p2.pointOnPlane(pt, 1E-6);
		if (side != side2) return true;
		//if (p1.pointOnPlane(pt, 1E-6) && p2.pointOnPlane(pt, 1E-6)) return true;
		return false;
	}

	bool isSilhouetteForRay(Ray& r) {
		return isSilhouetteForPos(r.origin);
	}

	//normal here points in the direction of the 'front' so edge must lie before vertex to be silhouette
	bool isSilhouetteForVertex(Vertex* v, bool &side, glm::dvec3 normal = glm::dvec3()) const {
		if (normal.length() != 0) {
			if (glm::dot(vertices[0]->pos, -normal) - glm::dot(v->pos, -normal) > -1E-7 &&
				glm::dot(vertices[1]->pos, -normal) - glm::dot(v->pos, -normal) > -1E-7) return false;
		}
		//if (v->id == vertices[0]->id || v->id == vertices[1]->id) return true;
		return isSilhouetteForPos(v->pos, side);
	}

	bool isInFrontOfEdge(Edge* e, Plane* plane) const {
		return  isInFrontOfVertex(e->vertices[0], plane) || isInFrontOfVertex(e->vertices[1], plane);
	}

	bool isInFrontOfVertex(Vertex* v, Plane* plane) const {
		return vertices[0]->isInFrontOfVertex(v, plane) || vertices[1]->isInFrontOfVertex(v, plane);
	}

	bool isSilhouetteForEdge(Edge* e) const {
		if (triangles.size() == 1) return true;
		Plane p1 = triangles[0]->getPlane();
		Plane p2 = triangles[1]->getPlane();
		int count = 0;
		for (Vertex* v : e->vertices) {
			bool side;
			if (isSilhouetteForPos(v->pos, p1, p2, side)) return true;
			if (side) count++;
		}
		if (count < 2) return true;
		return false;
	}

	bool isSilhouetteForPrim(std::vector<glm::dvec3>& pts) const {
		if (triangles.size() == 1) return true;
		Plane p1 = triangles[0]->getPlane();
		Plane p2 = triangles[1]->getPlane();
		int count = 0;
		for (glm::dvec3& pt : pts) {
			bool side;
			if (isSilhouetteForPos(pt, p1, p2, side)) return true;
			if (side) count++;
		}
		if (count < pts.size()) return true;
		return false;
	}

	bool convex() {
		if (triangles.size() == 2) {
			glm::dvec3 vNotOnEdge1pos;
			glm::dvec3 vNotOnEdge2pos;
			for (int i = 0; i < 3; i++) {
				int vNotOnEdge1 = triangles[0]->vertices[i]->id;
				int vNotOnEdge2 = triangles[1]->vertices[i]->id;
				if (vNotOnEdge1 != vertices[0]->id && vNotOnEdge1 != vertices[1]->id) vNotOnEdge1pos = triangles[0]->vertices[i]->pos;
				if (vNotOnEdge2 != vertices[0]->id && vNotOnEdge2 != vertices[1]->id) vNotOnEdge2pos = triangles[1]->vertices[i]->pos;
			}
			double concaveConvex = glm::dot(vNotOnEdge2pos - vNotOnEdge1pos, triangles[0]->normal);
			if (fabs(concaveConvex) < 1E-10 || concaveConvex > 0) return false; // in plane or concave
		}
		return true;
	}

	struct cmp_ptr
	{
		bool operator()(const Edge* lhs, const Edge* rhs) const
		{
			return lhs->id < rhs->id;
		}
	};

	uint64_t getKey() {
		return (uint64_t) vertices[0]->id << 32 | (uint64_t) vertices[1]->id;
	};

};

//struct EdgeComp
//{
//	bool operator()(const Edge* lhs, const Edge* rhs) const { 
//		return  lhs->vertices[0]->id < rhs->vertices[0]->id ||
//				(lhs->vertices[0]->id == rhs->vertices[0]->id &&
//				lhs->vertices[1]->id < rhs->vertices[1]->id);	}
//};

//struct EdgeSet {
//	int v0;
//	int v1;
//	mutable int id = -1;
//	bool operator<(const EdgeSet& e) const {
//		return  e.v0 < this->v0 || (e.v0 == this->v0 && e.v1 < this->v1);
//	}
//};
#endif