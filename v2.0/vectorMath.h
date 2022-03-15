#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>

static double dot_fd(glm::vec3 floatVec, glm::dvec3 doubleVec) {
	return glm::dot(glm::dvec3(floatVec), doubleVec);
};

static double dot_fd(glm::dvec3 doubleVec, glm::vec3 floatVec) {
	return glm::dot(glm::dvec3(floatVec), doubleVec);
};
