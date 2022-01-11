#pragma once
#include "camera.h"

struct Orthocamera : public Camera {

	float w, h;

	Orthocamera::Orthocamera(float w, float h, float far = 30.f, float near = 0.1f)
		: Camera{ near, far }, w(w), h(h) {};

	virtual Camera* makeCopy() override
	{
		Camera* ptr = new Orthocamera(*this);
		return ptr;
	};

	virtual glm::mat4 pMatrix() const override
	{
		return glm::ortho(-w, w, -h, h, near, far);
	};

	virtual Ray pixRayDirection(glm::vec2 pos) override
	{
		glm::vec3 camRay = glm::vec3(pos.x*w, pos.y*h, 0);
		glm::vec3 worldSpaceRayStart = invview * glm::vec4(camRay, 0.0f);
		//if (worldSpaceRayStart.length() > 0) worldSpaceRayStart = glm::normalize(worldSpaceRayStart);
		return Ray(position + worldSpaceRayStart, position + worldSpaceRayStart + forward);
	};
};