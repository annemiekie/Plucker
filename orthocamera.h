#pragma once
#include "camera.h"

struct Orthocamera : public Camera {

	float w, h;

	Orthocamera::Orthocamera(float w, float h, float near, float far)
		: Camera{ near, far }, w(w), h(h) {};

	virtual glm::mat4 pMatrix() const
	{
		return glm::ortho(-w, w, -h, h, near, far);
	};

	virtual Ray pixRayDirection(glm::vec2 pos)
	{
		glm::vec3 camRay = glm::vec3(pos.x*w, pos.y*h, 0);
		glm::vec3 worldSpaceRayStart = glm::normalize(invview * glm::vec4(camRay, 0.0f));
		return Ray(position + forward, position + worldSpaceRayStart);
	};
};