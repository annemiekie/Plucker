#pragma once
#include "plane.h"

class Polygone : public virtual Plane {
public:
	std::vector<Ray> edges;
	std::vector<glm::dvec3> vertices;

	Polygone(glm::dvec3 n, float c) : Plane(n, c) { };
	Polygone() {};
	~Polygone() {};

	Polygone(std::vector<glm::dvec3>& vertices, glm::dvec3 voppo) 
		: Plane(vertices, voppo), vertices(vertices) {
		makeDirectedEdges();
	}

	void makeDirectedEdges() {
		edges.resize(vertices.size());
		glm::dvec3 center;
		for (glm::dvec3& v : vertices) center += v;
		center /= (float)vertices.size();
		Ray centerInvNormal = Ray(normal + center, center);
		for (int i = 0; i < vertices.size(); i++) {
			edges[i] = Ray(vertices[i], vertices[(i + 1) % 4]);
			if (centerInvNormal.side(edges[i])) edges[i].inverseDir();
		}
	}

	glm::dvec3 getMin() {
		glm::dvec3 minim = glm::dvec3(INFINITY);
		for (glm::dvec3& v : vertices) {
			for (int i = 0; i < 3; i++) {
				if (v[i] < minim[i]) minim[i] = v[i];
				//if (v.x <= minim.x && v.y <= minim.y && v.z <= minim.z) minim = v;
			}
		}
		return minim;
	}

	glm::dvec3 getMax() {
		glm::dvec3 maxim = glm::vec3(-INFINITY);
		for (glm::dvec3& v : vertices) {
			for (int i = 0; i < 3; i++) {
				if (v[i] > maxim[i]) maxim[i] = v[i];
			}
		}
		return maxim;
	}

	bool inBounds(const Ray& r, float thres = 0) {
		int count = 0;
		for (Ray& e : edges) {
			double x = e.sideVal(r);
			if (x < -thres) return false;
			if (fabs(x) < thres) count++;
		}
		if (count == 4) return false;
		return true;
	}

	bool inBounds(const glm::dvec3& pt, float thres = 0) {
		Ray r(pt + normal, pt);
		return (inBounds(r, thres));
	}

	bool intersectsPlaneFromLines(std::vector<Ray>& lines) {
		if (lines.size() == 1) {
			Ray r = lines[0];
			if (inBounds(r)) return true;
			r.inverseDir();
			if (inBounds(r)) return true;
			return false;
		}
		glm::dvec3 pt1 = rayIntersection(lines[0]);
		glm::dvec3 pt2 = rayIntersection(lines[1]);
		Ray r12(pt1, pt2);
		for (Ray& e : edges) {
			glm::dvec3 intersect = e.pointOfintersectWithRay(r12);
			if (inBounds(intersect)) return true;
		}
		return false;
	}

};