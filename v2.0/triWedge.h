#pragma once
#include "plane.h"

class TriWedge : public Plane {
public:
	glm::vec3 vec1;
	glm::vec3 vec2;
	glm::vec3 vMain;

    TriWedge(glm::vec3 n, float c) : Plane(n,c) { };

	TriWedge(glm::vec3 v1main, glm::vec3 v2, glm::vec3 v3, glm::vec3 voppo) : Plane(v1main, v2, v3, voppo), vMain(v1main) {
		vec1 = v2 - v1main;
		vec2 = v3 - v1main;
	};

	bool inBounds(glm::vec3& pt) {	
		return glm::dot(glm::cross(vec1, pt - vMain), glm::cross(vec2, pt - vMain)) <= 0;
		//bool dot1 = glm::dot(vec1, pt - vMain) > 0;
		//bool dot2 = glm::dot(vec2, pt - vMain) > 0;
		//return dot1 == dot2;
	}

	virtual bool raySegmentIntersection(glm::vec3 v1, glm::vec3 v2, glm::vec3& pt) override {
		if (!Plane::raySegmentIntersection(v1, v2, pt)) return false;
		return inBounds(pt);
	}
};