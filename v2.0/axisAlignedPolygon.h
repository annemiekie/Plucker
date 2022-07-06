#pragma once
#include "polygon.h"
#include "axisAlignedPlane.h"
#include "segment2D.h"

class AxisAlignedPolygon : virtual public Polygone, public AxisAlignedPlane {
public:
	std::vector<glm::dvec2> vertices2D;

	AxisAlignedPolygon() {};
	~AxisAlignedPolygon() {};
	AxisAlignedPolygon(AxisAlignedPolygon* polygon) : Plane(polygon), AxisAlignedPlane(polygon) {
		vertices = polygon->vertices;
		s = polygon->s;
		edges = polygon->edges;
		vertexEdges = polygon->vertexEdges;
		vertices2D = polygon->vertices2D;
	};

	AxisAlignedPolygon(AxisAlignedPlane* plane) : Plane(plane), AxisAlignedPlane(plane) {};

	AxisAlignedPolygon(glm::dvec3 point, glm::dvec3 normal) : Plane(point, normal), AxisAlignedPlane(point, normal) {	};

	AxisAlignedPolygon(std::vector<glm::dvec3> vertices, glm::ivec3 normal) : Plane(vertices, normal), AxisAlignedPlane(vertices[0], normal), Polygone(vertices, normal) {		
		for (int i = 0; i < vertices.size(); i++) {
			vertices2D[i] = filter3Dto2D(vertices[i]);
		}
	}

	glm::dvec2 filter3Dto2D(glm::dvec3 vec3) {
		glm::dvec2 vec2;
		int c2 = 0;
		for (int c3 = 0; c3 < 3; c3++) {
			if (c3 != constantCoord) {
				vec2[c2] = vec3[c3];
				c2++;
			}
		}
		return vec2;
	}

	bool intersectionInBounds(Segment2D& rseg) {
		for (int i = 0; i < vertices2D.size(); i++) {
			glm::dvec2 v1 = vertices2D[i];
			glm::dvec2 v2 = vertices2D[(i+1) % vertices2D.size()];
			Segment2D vseg = { v1, v2 };
			if (vseg.intersectSegmSegm(rseg)) return true;
		}
		return false;
	}

	bool intersectsPlaneFromLines(std::vector<Ray>& lines) {
		for (Ray& r : lines) if (lineInBounds(r)) return true;
		glm::dvec2 pt1 = filter3Dto2D(rayIntersection(lines[0]));
		glm::dvec2 pt2 = filter3Dto2D(rayIntersection(lines[1]));
		Segment2D pts = { pt1, pt2 };
		return intersectionInBounds(pts);
	}

	double getMaxSize() {
		glm::vec2 size = getMaxAA() - getMinAA();
		return std::max(size[0], size[1]);
	}

	glm::dvec2 getMinAA() {
		glm::dvec2 minim = glm::dvec2(INFINITY);
		for (glm::dvec2& v : vertices2D) {
			for (int i = 0; i < 2; i++) {
				if (v[i] < minim[i]) minim[i] = v[i];
				//if (v.x <= minim.x && v.y <= minim.y && v.z <= minim.z) minim = v;
			}
		}
		return minim;
	}

	glm::dvec2 getMaxAA() {
		glm::dvec2 maxim = glm::vec2(-INFINITY);
		for (glm::dvec2& v : vertices2D) {
			for (int i = 0; i < 2; i++) {
				if (v[i] > maxim[i]) maxim[i] = v[i];
			}
		}
		return maxim;
	}


};