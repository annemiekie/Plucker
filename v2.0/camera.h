#ifndef CAMERA_H
#define CAMERA_H

// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

// Library for window creation and event handling
#include <GLFW/glfw3.h>
#include "ray.h"

struct Camera
{

	glm::vec3 position;
	glm::mat4 invview;
	glm::vec3 forward;
	glm::vec3 up;

	float near;
	float far;


	Camera(float near = 0.1f, float far = 30.f)
		: position(glm::vec3(0, 0, 0))
		, forward(glm::vec3(0, 0, -1))
		, up(glm::vec3(0, 1, 0))
		, near(near)
		, far(far)
	{
		invview = invVMatrix();
	};

	virtual Camera* makeCopy() {
		Camera* ptr = new Camera(*this);
		return ptr;
	};

	virtual Ray pixRayDirection(glm::vec2 pos) {
		return Ray();
	};

	void setPositionAndForward(glm::vec3 pos, glm::vec3 lookat) {
		position = pos;
		forward = lookat - pos;
		invview = invVMatrix();
	}

	void setPosition(glm::vec3 pos) {
		position = pos;
		invview = invVMatrix();
	}

	void setForward(glm::vec3 fw) {
		forward = fw;
		invview = invVMatrix();
	}

	virtual glm::mat4 invVMatrix() const
	{
		return glm::inverse(vMatrix());
	}

	virtual glm::mat4 vMatrix() const
	{
		return glm::lookAt(position, position + forward, up);
	}

	virtual glm::mat4 pMatrix() const
	{
		return glm::mat4(0);
	};

	virtual glm::mat4 vpMatrix() const
	{
		return pMatrix() * vMatrix();
	}

};

#endif