#pragma once
// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

// Library for window creation and event handling
#include <GLFW/glfw3.h>

struct PreviewCamera {

	float y_rot = 0.f;
	float x_rot = 0.f;

	glm::mat4 world = glm::mat4(1);
	float radius = 6.f;
	glm::mat4 view = glm::translate(glm::mat4(1), glm::vec3(0.f, 0.f, -radius));

	glm::vec3 up;

	float near;
	float far;

	float fov;
	float aspect;
	float halfscreen;

	glm::vec3 position;
	glm::mat4 invview;
	glm::vec3 lookat;


	PreviewCamera::PreviewCamera()
		: up(glm::vec3(0, 1, 0))
		, lookat(glm::vec3(0, 0, 0))
		, fov(glm::pi<float>() / 4.f)
		, aspect(1.f)
		, near(0.1f)
		, far(30.f)
		, halfscreen(std::tan(fov / 2.f))
	{
		invview = invVMatrix();
		position = radius * glm::vec3(1, 0, 0);
	}

	Ray pixRayDirection(glm::vec2 pos)
	{
		glm::vec3 camRay = glm::normalize(glm::vec3(pos.x * halfscreen, pos.y * halfscreen * aspect, -1));
		glm::mat4 inv = glm::inverse(view*world);
		glm::vec3 worldSpaceRayDirection = glm::normalize(inv * glm::vec4(camRay, 0.0f));
		glm::vec3 posi = inv[3];
		return Ray(posi + worldSpaceRayDirection, posi);
	};

	glm::mat4 pMatrix() const
	{
		return glm::perspective(fov, aspect, near, far);
	}

	glm::mat4 invVMatrix() const
	{
		return glm::inverse(vMatrix());
	}

	glm::mat4 vMatrix() const
	{
		return glm::lookAt(position, lookat, up);
	}

	glm::mat4 vpMatrix() const
	{
		return pMatrix() * view * world;
		//return pMatrix() * vMatrix();
	}

	glm::mat4 vpmMatrix(glm::vec3 translate) const
	{
		glm::mat4 model = glm::translate(glm::mat4(1), translate);
		return pMatrix() * view * world * model;
		//return pMatrix() * vMatrix() * model;

	}

	void updateRadius(float radiusChange) {
		if (radiusChange == 0) return;
		else {
			radiusChange /= 4.f;
			radius += radiusChange;
			if (radius <= 0.f) radius = 6.f;
			glm::mat4 camToWorld = glm::inverse(view);

			glm::vec3 viewpos = camToWorld[3];
			glm::vec3 norm = glm::normalize(viewpos);
			view = glm::translate(view, norm * radiusChange);
			
		}
	}

	void updateDirection(glm::vec2 rotSpeed, int width, int height)
	{
		// Tweak these values to change the sensitivity
		float scale_x = rotSpeed.x / (float)width;
		float scale_y = rotSpeed.y / (float)height;
		float rotSpeedTot = 350.0f;
		float y_rot_clamp = 89.9999f;
		float rotsigny = rotSpeed.y < 0 ? 1 : -1;


		float rot = rotSpeedTot * scale_x;
		world = glm::rotate(world, (float)glm::radians(rot), up);
		glm::mat4 rotationM = glm::rotate(glm::mat4(), (float)glm::radians(-rot), up);
		position = (glm::mat3)rotationM * position;
		x_rot += rot;

		rot = rotSpeedTot * -scale_y;
		if (abs(y_rot + rot) > y_rot_clamp) rot = rotsigny *y_rot_clamp - y_rot;

		rotationM = glm::rotate(glm::mat4(), (float)glm::radians(rot), glm::normalize(glm::vec3(-position.y, position.x, 0)));
		position = (glm::mat3)rotationM * position;

		view = rotate(view, (float)glm::radians(rot), glm::vec3(1, 0, 0));
		y_rot += rot;
	}
};

float radiusChange = 0;

glm::vec3 camSpeed = glm::vec3(0, 0, 0);
glm::vec2 camRotSpeed = glm::vec2(0, 0);

bool      rightMousePressed = false;
glm::vec2 cursorPos = glm::vec2(0.0, 0.0);



// Mouse button handle function
void camMouseButtonHandler(int button, int action)
{
	if (button == GLFW_MOUSE_BUTTON_1 && action == GLFW_PRESS)   rightMousePressed = true;
	if (button == GLFW_MOUSE_BUTTON_1 && action == GLFW_RELEASE) rightMousePressed = false;
}

void camCursorPosHandler(double xpos, double ypos, int width, int height)
{

	if (!rightMousePressed) {
		camRotSpeed = glm::vec2(0, 0);
	}
	else {
		camRotSpeed.x = xpos - cursorPos.x;
		camRotSpeed.y = cursorPos.y - ypos;
	}

	cursorPos.x = (float)xpos;
	cursorPos.y = (float)ypos;
}

void camScrollHandler(double yoffset) {
	radiusChange = yoffset;
}

void updateCamera(PreviewCamera& camera, int width, int height)
{
	//camera.updatePosition(camSpeed);
	camera.updateDirection(camRotSpeed, width, height);
	camera.updateRadius(radiusChange);
	radiusChange = 0;
	camRotSpeed = glm::vec2(0, 0);
}

