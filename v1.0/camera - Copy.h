// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "ray.h"

// Library for window creation and event handling
#include <GLFW/glfw3.h>

class Camera
{	
public:
	glm::vec3 position;
	glm::vec3 forward;
	glm::vec3 up;

	glm::mat4 world = glm::mat4(1);
	float radius = 3.f;
	glm::mat4 view = glm::translate(glm::mat4(1), glm::vec3(0, 0, -radius));
	
	float near;
	float far;

	float     fov;
	float     aspect;

	float y_rot = 0.f;
	float x_rot = 0.f;
	float halfscreen = 0.f;
	glm::mat4 invVMatrix;

	Camera::Camera()
		: position(glm::vec3(0, 0, 0))
		, forward(glm::vec3(0, 0, -1))
		, up(glm::vec3(0, 1, 0))
		, fov(glm::pi<float>() / 4.f)
		, aspect(1.f)
		, near(0.1f)
		, far(30.f)
	{
		invVMatrix = glm::inverse(vMatrix());
		halfscreen = std::tan(fov / 2.f);
	};

	glm::mat4 vMatrix() const
	{
		return glm::lookAt(position, position + forward, up);
	};

	glm::mat4 pMatrix() const
	{
		return glm::perspective(fov, aspect, near, far);
	};

	glm::mat4 oMatrix(float w, float h) const
	{
		return glm::ortho(-2.f, 2.f, -2.f, 2.f, near, far);
	};

	glm::mat4 vpMatrix() const
	{
		return pMatrix() * vMatrix();
	};

	glm::mat4 vpmMatrix() const
	{
		glm::mat4 model = glm::translate(glm::mat4(1), glm::vec3(0.0f, -0.5f, 0.0f));
		return pMatrix() * view * world * model;
	};

	glm::mat4 vomMatrix(int w, int h) const
	{
		glm::mat4 model = glm::translate(glm::mat4(1), glm::vec3(0.0f, -0.5f, 0.0f));
		return oMatrix(1.f * w, 1.f * h) * view * world * model;
	};

	void updatePosition(glm::vec3 speed)
	{
		forward = -position;
		forward = glm::normalize(forward);
		up = glm::normalize(up);
		glm::vec3 right = glm::normalize(glm::cross(forward, up));

		position += speed.z * forward + speed.y * up + speed.x * right;
	};

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
	};
	
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
	};
};

float radiusChange = 0;

glm::vec3 camSpeed    = glm::vec3(0, 0, 0);
glm::vec2 camRotSpeed = glm::vec2(0, 0);

bool      rightMousePressed = false;
glm::vec2 cursorPos = glm::vec2(0.0, 0.0);


// Key handle function
void cameraKeyboardHandler(int key, int action)
{
	const float sp = 0.05f;

	switch (key) 
	{
	case GLFW_KEY_A:
		if (action == GLFW_PRESS)   camSpeed.x =  -sp;
		if (action == GLFW_RELEASE) camSpeed.x =  0.0;
		break;
	case GLFW_KEY_D:
		if (action == GLFW_PRESS)   camSpeed.x =  sp;
		if (action == GLFW_RELEASE) camSpeed.x =  0.0;
		break;
	case GLFW_KEY_W:
		if (action == GLFW_PRESS)   camSpeed.z =  sp;
		if (action == GLFW_RELEASE) camSpeed.z =  0.0;
		break;
	case GLFW_KEY_S:
		if (action == GLFW_PRESS)   camSpeed.z = -sp;
		if (action == GLFW_RELEASE) camSpeed.z =  0.0;
		break;
	case GLFW_KEY_R:
		if (action == GLFW_PRESS)   camSpeed.y =  sp;
		if (action == GLFW_RELEASE) camSpeed.y =  0.0;
		break;
	case GLFW_KEY_F:
		if (action == GLFW_PRESS)   camSpeed.y = -sp;
		if (action == GLFW_RELEASE) camSpeed.y =  0.0;
		break;
	default:
		break;
	}
}


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
	} else {
	//	camRotSpeed = glm::vec2(xpos, ypos) - cursorPos;

		// Calculate the distance, the mouse was moved,
		// between the last and the current frame
		camRotSpeed.x = xpos - cursorPos.x;
		camRotSpeed.y = cursorPos.y - ypos;
	}

	cursorPos.x = (float)xpos;
	cursorPos.y = (float)ypos;
}

void camScrollHandler(double yoffset) {
	radiusChange = yoffset;
}

void updateCamera(Camera& camera, int width, int height) 
{
	//camera.updatePosition(camSpeed);
	camera.updateDirection(camRotSpeed, width, height);
	camera.updateRadius(radiusChange);
	radiusChange = 0;
	camRotSpeed = glm::vec2(0, 0);
}

