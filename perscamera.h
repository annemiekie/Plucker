#pragma once
#include "camera.h"

struct PersCamera : public Camera {

	float     fov;
	float     aspect;
	float halfscreen = 0.f;

	PersCamera(float near, float far)
		: Camera{near, far}, fov(glm::pi<float>() / 4.f), aspect(1.f)
	{
		halfscreen = std::tan(fov / 2.f);
	}

	virtual Ray pixRayDirection(glm::vec2 pos)
	{
		glm::vec3 camRay = glm::normalize(glm::vec3(pos.x * halfscreen, pos.y * halfscreen * aspect, -1));
		glm::vec3 worldSpaceRayDirection = glm::normalize(invview * glm::vec4(camRay, 0.0f));
		return Ray(position + worldSpaceRayDirection, position);
	};

	virtual glm::mat4 pMatrix() const
	{
		return glm::perspective(fov, aspect, near, far);
	}

}; 
