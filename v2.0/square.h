#pragma once
#include "plane.h"

class Square : public Plane {
public:
	std::vector<Ray> quadLines = std::vector<Ray>(4);
	std::vector<glm::dvec3> cornerPoints = std::vector<glm::dvec3>(4);

	Square(glm::dvec3 n, float c) : Plane(n, c) { };
	Square() {};
	~Square() {};

	Square(glm::dvec3 v1main, glm::dvec3 v2, glm::dvec3 v3, glm::dvec3 voppo) : Plane(v1main, v2, v3, voppo) {
		cornerPoints[0] = v1main;
		cornerPoints[1] = v2;
		cornerPoints[2] = v2 + v3 - v1main;
		cornerPoints[3] = v3;

		glm::dvec3 center;
		for (glm::dvec3& c : cornerPoints) center += c;
		center /= 4.;
		Ray centerInvNormal = Ray(normal + center, center);

		for (int i = 0; i < 4; i++) {
			quadLines[i] = Ray(cornerPoints[i], cornerPoints[(i + 1) % 4]);
			if (centerInvNormal.side(quadLines[i])) quadLines[i].inverseDir();
			//float x = centerInvNormal.sideVal(quadLines[i]);
		}
	}

	glm::dvec3 getMin() {
		glm::dvec3 minim = glm::dvec3(INFINITY);
		for (glm::dvec3& cp : cornerPoints) {
			if (cp.x <= minim.x && cp.y <= minim.y && cp.z <= minim.z) minim = cp;
		}
		return minim;
	}

	glm::dvec3 getMax() {
		glm::dvec3 maxim = glm::vec3(-INFINITY);
		for (glm::dvec3& cp : cornerPoints) {
			if (cp.x >= maxim.x && cp.y >= maxim.y && cp.z >= maxim.z) maxim = cp;
		}
		return maxim;
	}

	bool inBounds(const Ray& r, float thres = 0) {
		int count = 0;
		for (Ray& ql : quadLines) {
			double x = ql.sideVal(r);
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
		//for (Ray& r : lines) {
		//	if (inBounds(r)) return true;
		//	r.inverseDir();
		//	if (inBounds(r)) return true;
		//}
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
		for (Ray& qr : quadLines) {
			glm::dvec3 intersect = qr.pointOfintersectWithRay(r12);
			if (inBounds(intersect)) return true;
		}
		return false;
	}

}; 

#pragma once
