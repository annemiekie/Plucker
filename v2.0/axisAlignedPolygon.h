#pragma once
#include "polygon.h"
#include "axisAlignedPlane.h"

class AxisAlignedPolygon : virtual public Polygone, public AxisAlignedPlane {
public:
	std::vector<glm::dvec2> vertices2D;

	AxisAlignedPolygon() {};
	~AxisAlignedPolygon() {};
	AxisAlignedPolygon(AxisAlignedPolygon* polygon) : Plane(polygon), AxisAlignedPlane(polygon) {
		vertices = polygon->vertices;
		edges = polygon->edges;
		vertices2D = polygon->vertices2D;
	};

	AxisAlignedPolygon(AxisAlignedPlane* plane) : Plane(plane), AxisAlignedPlane(plane) {};


	AxisAlignedPolygon(glm::dvec3 point, glm::dvec3 normal) : Plane(point, normal), AxisAlignedPlane(point, normal) {	};

	AxisAlignedPolygon(std::vector<glm::dvec3> vertices, glm::ivec3 normal) : Plane(vertices, normal), AxisAlignedPlane(vertices[0], normal), Polygone(vertices, normal) {
		
		for (int i = 0; i < vertices.size(); i++) {
			int c2 = 0;
			for (int c3 = 0; c3 < 3; c3++) {
				if (c3 != constantCoord) {
					vertices2D[i][c2] = vertices[i][c3];
					c2++;
				}
				
			}
		}
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