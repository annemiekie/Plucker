#pragma once
#include "Ray.h"

class Line4 : public Ray {
public:
	int index = -1;
	bool checkboth = true;
	std::set<int> prims;
	bool inplane = false;
	Line4() : Ray() {};
	Line4(const Ray& ray) : Ray(ray) {};
	Line4(glm::dvec3 from, glm::dvec3 to) : Ray(from, to) {};
	Line4(Eigen::VectorXd plucker, bool reverse) : Ray(plucker, reverse) {	};
};