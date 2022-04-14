#pragma once
#include "square.h"
#include "axisAlignedPolygon.h"
#include "axisAlignedPlane.h"

class AxisAlignedSquare : public Square, public AxisAlignedPolygon {
public:

	//AxisAlignedPolygon(std::vector<glm::dvec3> vertices, glm::ivec3 normal) : AxisAlignedPolygon(vertices[0], normal) {

	//}
	glm::dvec2 minm = { INFINITY, INFINITY };
	glm::dvec2 maxm = { -INFINITY, -INFINITY };

	AxisAlignedSquare(std::vector<glm::dvec3>& vertices, glm::ivec3 normal) : Plane(vertices[0], normal), AxisAlignedPolygon(vertices, normal) {
		getMinMaxFromPoints(vertices);
	}

	AxisAlignedSquare(glm::dvec2& minm, glm::dvec2& maxm, glm::dvec3 point, glm::ivec3 normal) 
		: Plane(point, normal), AxisAlignedPolygon(point, normal), minm(minm), maxm(maxm) {
		makeVerticesFromMinMax();
		makeDirectedEdges();

	}

	AxisAlignedSquare(std::vector<Ray>& rays, AxisAlignedPlane* plane) : Plane(plane), AxisAlignedPolygon(plane) {
		std::vector<glm::dvec3> points;
		for (Ray& r : rays) points.push_back(rayIntersection(r));
		getMinMaxFromPoints(points);
		makeVerticesFromMinMax();
		makeDirectedEdges();
	};

	AxisAlignedSquare(std::vector<glm::dvec3>& points, AxisAlignedPlane* plane) : Plane(plane), AxisAlignedPolygon(plane) {
		getMinMaxFromPoints(points);
		makeVerticesFromMinMax();
		makeDirectedEdges();
	};

	AxisAlignedSquare(AxisAlignedPolygon* polygon) : Plane(polygon), AxisAlignedPolygon(polygon) {
		getMinMaxFromPoints(polygon->vertices);
		makeVerticesFromMinMax();
		makeDirectedEdges();
	}

	void getMinMaxFromPoints(std::vector<glm::dvec3>& points) {
		for (glm::dvec3 v : points) {
			int c2 = 0;
			for (int c3 = 0; c3 < 3; c3++) {
				if (c3 != constantCoord) {
					if (v[c3] < minm[c2]) minm[c2] = v[c3];
					if (v[c3] > maxm[c2]) maxm[c2] = v[c3];
					c2++;
				}
			}
		}
	}

	void makeVerticesFromMinMax() {
		vertices2D = { {minm[0], minm[1] }, {maxm[0], minm[1]}, {maxm[0], maxm[1]}, {minm[0], maxm[1]} };
		vertices.resize(vertices2D.size());
		for (int i = 0; i < vertices2D.size(); i++) {
			int c2 = 0;
			vertices[i] = aaPosition;
			for (int c3 = 0; c3 < 3; c3++) {
				if (c3 != constantCoord) {
					vertices[i][c3] = vertices2D[i][c2];
					c2++;
				}
			}
		}

	};

	bool noIntersection(AxisAlignedSquare& other) {
		for (int i = 0; i < 2; i++) {
			if (minm[i] > other.maxm[i]) return true;
			if (maxm[i] < other.minm[i]) return true;
		}
		return false;
	};

};