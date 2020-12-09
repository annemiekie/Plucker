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

	float     fov;
	float     aspect;

	PreviewCamera::PreviewCamera()
		: up(glm::vec3(0, 1, 0))
		, fov(glm::pi<float>() / 4.f)
		, aspect(1.f)
		, near(0.1f)
		, far(30.f)
	{}

	glm::mat4 pMatrix() const
	{
		return glm::perspective(fov, aspect, near, far);
	}

	glm::mat4 vpmMatrix() const
	{
		glm::mat4 model = glm::translate(glm::mat4(1), glm::vec3(0.f, -1.f, 0.f));
		return pMatrix() * view * world *model;
	}

	glm::mat4 vpmMatrix(glm::vec3 translate) const
	{
		glm::mat4 model = glm::translate(glm::mat4(1), glm::vec3(0.f, -1.f, 0.f));
		model = glm::translate(model, translate);
		return pMatrix() * view * world * model;
	}

	void updateRadius(float radiusChange) {
		if (radiusChange == 0) return;
		else {
			radiusChange /= 4.f;
			radius += radiusChange;
			glm::mat4 camToWorld = glm::inverse(view);

			glm::vec3 viewpos = camToWorld[3];
			glm::vec3 norm = glm::normalize(viewpos);
			view = glm::translate(view, norm * radiusChange);
		}
	}

	void updateDirection(glm::vec2 rotSpeed, int width, int height)
	{
		// Tweak these values to change the sensitivity
		float scale_x = abs(rotSpeed.x) / (float)width;
		float scale_y = abs(rotSpeed.y) / (float)height;
		float rotSpeedTot = 350.0f;
		float y_rot_clamp = 89.9999f;

		// Horizontal rotation (on the Y-axis)
		// This is simple because no clamping is needed
		if (rotSpeed.x < 0)
		{
			// As discussed earlier, the entire world is rotated
			world = glm::rotate(world, (float)glm::radians(-rotSpeedTot * scale_x), up);
			x_rot -= rotSpeedTot * scale_x;
		}
		else if (rotSpeed.x > 0)
		{
			world = glm::rotate(world, (float)glm::radians(rotSpeedTot * scale_x), up);
			x_rot += rotSpeedTot * scale_x;
		}

		// The user wants to rotate the camera this much
		float rot = rotSpeedTot * scale_y;

		if (rotSpeed.y < 0)
		{
			// Upper rotation limit (+90 deg)
			if (y_rot + rot > y_rot_clamp)
				rot = y_rot_clamp - y_rot;

			view = rotate(view, (float)glm::radians(rot), glm::vec3(1, 0, 0));
			y_rot += rot;
		}
		else if (rotSpeed.y > 0)
		{
			// Limit the rotation in the other direction too (-90 deg)
			if (y_rot - rot < -y_rot_clamp)
				rot = y_rot + y_rot_clamp;

			view = glm::rotate(view, (float)glm::radians(-rot), glm::vec3(1, 0, 0));
			y_rot -= rot;
		}
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

