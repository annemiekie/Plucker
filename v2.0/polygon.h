#pragma once
#include "plane.h"

class Polygone : public virtual Plane {
public:
	std::vector<Ray> edges;
	std::vector<glm::dvec3> vertices;
	std::vector<std::vector<glm::dvec3>> vertexEdges;
	int s = 0;
	GLuint vao;

	Polygone(glm::dvec3 n, float c) : Plane(n, c) { };
	Polygone() {};
	~Polygone() {};

	Polygone(std::vector<glm::dvec3>& vertices, glm::dvec3 voppo) 
		: Plane(vertices, voppo), vertices(vertices) {
		s = vertices.size();
		makeDirectedEdges();
		makeVertexEdges();
	}

	void makeVertexEdges() {
		vertexEdges.resize(vertices.size());
		for (int i = 0; i < vertices.size(); i++) {
			vertexEdges[i] = { vertices[(i - 1 + vertices.size()) % vertices.size()], vertices[i] };
		}
	}

	void makeDirectedEdges() {
		edges.resize(vertices.size());
		glm::dvec3 center;
		for (glm::dvec3& v : vertices) center += v;
		center /= (float)vertices.size();
		Ray centerInvNormal = Ray(normal + center, center);
		for (int i = 0; i < vertices.size(); i++) {
			edges[i] = Ray(vertices[i], vertices[(i + 1) % s]);
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

	bool inBounds(const Ray& r, double thres = 1E-8) {
		int count = 0;
		for (Ray& e : edges) {
			double x = e.sideVal(r);
			if (x < -thres) return false;
			if (fabs(x) < thres) count++;
		}
		if (count == 4) return false;
		return true;
	}

	bool inBounds(const glm::dvec3& pt, double thres = 1E-8) {
		Ray r(pt + normal, pt);
		return (inBounds(r, thres));
	}

	bool lineInBounds(Ray& r, double thres = 1E-8) {
		if (inBounds(r)) return true;
		r.inverseDir();
		if (inBounds(r)) return true;
		return false;
	}

	bool planeIntersectionFromLines(std::vector<Ray>& lines) {
		for (Ray& r : lines) if (lineInBounds(r)) return true;
		glm::dvec3 pt1 = rayIntersection(lines[0]);
		glm::dvec3 pt2 = rayIntersection(lines[1]);
		Ray r12(pt1, pt2);
		for (Ray& e : edges) {
			double depth;
			glm::dvec3 intersect = e.pointOfintersectWithRay(r12, depth);
			if (inBounds(intersect)) return true;
		}
		return false;
	}

	void vaoGeneration() {
		std::vector<glm::dvec3> lineSegments;
		for (int i = 0; i < vertices.size(); i++) {
			lineSegments.push_back(vertices[i]);
			lineSegments.push_back(vertices[(i+1)%vertices.size()]);
		}
		GLuint lineVBO;
		glGenBuffers(1, &lineVBO);
		glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
		glBufferData(GL_ARRAY_BUFFER, lineSegments.size() * sizeof(glm::dvec3), &lineSegments[0].x, GL_STATIC_DRAW);

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(glm::dvec3), (void*)0);
	}

	void draw() {
		glBindVertexArray(vao);
		glDrawArrays(GL_LINES, 0, vertices.size());
	};


};