#pragma once
#include "Ray.h"

class Line4 : public Ray {
public:
	int index = -1;
	bool checkboth = true;
	Line4() : Ray() {};
	Line4(glm::vec3 from, glm::vec3 to) : Ray(from, to) {};
	Line4(Eigen::VectorXd plucker, bool reverse) : Ray(plucker, reverse) {	};
};