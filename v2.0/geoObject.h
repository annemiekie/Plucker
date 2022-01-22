#pragma once

#ifndef GEOOBJECT_H
#define GEOOBJECT_H

#include <GL/glew.h>

#include "ray.h"
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>

class GeoObject {
public:
	GLuint vao;
	glm::vec3 center = glm::vec3(0);
	GeoObject() {};

	GeoObject(glm::vec3 center) : center(center) {};
	
	virtual bool intersect(const Ray& ray, glm::vec3& start, glm::vec3& end, bool getcolor = false, glm::vec3& maindir = glm::vec3(0), glm::vec3& color = glm::vec3()) {
		return false;
	};

	virtual void draw() {};

};

#endif