#pragma once
#include <glm/gtc/quaternion.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

class Window;

class Trackball {
public:
	// NOTE(Mathijs): field of view in radians! (use glm::radians(...) to convert from degrees to radians).
	Trackball(Window* pWindow, float fovy, float distanceFromLookAt = 4.0f, float rotationX = 0.0f, float rotationY = 0.0f);
	Trackball(Window* pWindow, float fovy, const glm::vec3& lookAt, float distanceFromLookAt = 4.0f, float rotationX = 0.0f, float rotationY = 0.0f);
	~Trackball() = default;

	static void printHelp();

	void disableTranslation();

	[[nodiscard]] glm::vec3 left() const;
	[[nodiscard]] glm::vec3 up() const;
	[[nodiscard]] glm::vec3 forward() const;

	[[nodiscard]] glm::vec3 position() const; // Position of the camera.
	[[nodiscard]] glm::vec3 lookAt() const; // Point that the camera is looking at / rotating around.
	[[nodiscard]] glm::mat4 viewMatrix() const;
	[[nodiscard]] glm::mat4 projectionMatrix() const;

	void setCamera(const glm::vec3 lookAt, const glm::vec3 rotations, const float dist); // Set the position and orientation of the camera.

private:

	void mouseButtonCallback(int button, int action, int mods);
	void mouseMoveCallback(const glm::vec2& pos);
	void mouseScrollCallback(const glm::vec2& offset);

private:
	const Window* m_pWindow;
	float m_fovy;
	bool m_canTranslate { true };

	glm::vec3 m_lookAt{ 0.0f }; // Point that the camera is looking at / rotating around.
	float m_distanceFromLookAt;
	glm::vec3 m_rotationEulerAngles{ 0 }; // Rotation as euler angles (in radians).

	glm::vec2 m_prevCursorPos; // Cursor position.
};