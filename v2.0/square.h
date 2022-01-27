#pragma once
#include "plane.h"

class Square : public Plane {
public:
	std::vector<Ray> quadLines = std::vector<Ray>(4);
	std::vector<glm::vec3> cornerPoints = std::vector<glm::vec3>(4);

	Square(glm::vec3 n, float c) : Plane(n, c) { };

	~Square() {};

	Square(glm::vec3 v1main, glm::vec3 v2, glm::vec3 v3, glm::vec3 voppo) : Plane(v1main, v2, v3, voppo) {
		cornerPoints[0] = v1main;
		cornerPoints[1] = v2;
		cornerPoints[2] = v2 + v3 - v1main;
		cornerPoints[3] = v3;

		glm::vec3 center;
		for (glm::vec3& c : cornerPoints) center += c;
		center /= 4.f;
		Ray centerInvNormal = Ray(normal + center, center);

		for (int i = 0; i < 4; i++) {
			quadLines[i] = Ray(cornerPoints[i], cornerPoints[(i + 1) % 4]);
			if (centerInvNormal.side(quadLines[i])) quadLines[i].inverseDir();
			//float x = centerInvNormal.sideVal(quadLines[i]);
		}
	}

	bool inBounds(const Ray& r, float thres = 0) {
		for (Ray& ql : quadLines) {
			float x = ql.sideVal(r);
			if (ql.sideVal(r) < -thres) return false;
		}
		return true;
	}

}; 

#pragma once
