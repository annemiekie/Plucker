// Library for OpenGL function loading
// Must be included before GLFW
#define GLEW_STATIC
#include <GL/glew.h>

// Library for window creation and event handling
#include <GLFW/glfw3.h>

// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>

// Library for loading an image
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

// Library for loading .OBJ model
#define TINYOBJLOADER_IMPLEMENTATION
//#include <tiny_obj_loader.h>

// Header for camera structure/functions
#include "perscamera.h"
#include "previewcamera.h"
#include "orthocamera.h"

#include "ray.h"
#include "vertex.h"
#include "rst.h"
#include "cube.h"
#include "models.h"
#include "shader.h"
#include "buffers.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <CImg.h>


// Configuration
const int width = 800;
const int height = 800;
int leafNumber = -1;

// Key handle function
void keyboardHandler(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	//cameraKeyboardHandler(key, action);
	if (action == GLFW_PRESS) {

		switch (key)
		{
		case GLFW_KEY_1:
			leafNumber++;
			break;
		case GLFW_KEY_2:
			break;
		default:
			break;
		}
	}
}

// Mouse button handle function
void mouseButtonHandler(GLFWwindow* window, int button, int action, int mods)
{
	camMouseButtonHandler(button, action);
}

void cursorPosHandler(GLFWwindow* window, double xpos, double ypos)
{
	camCursorPosHandler(xpos, ypos, width, height);
}

void scrollHandler(GLFWwindow* window, double xoffset, double yoffset) {
	camScrollHandler(yoffset);
}

// OpenGL debug callback
void APIENTRY debugCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {
	if (severity != GL_DEBUG_SEVERITY_NOTIFICATION) {
		std::cerr << "OpenGL: " << message << std::endl;
	}
}

void makeSample(glm::ivec2 resolution, Camera& cam, GLfloat* pixels, int triNum, RaySpaceTree* rst) {

	glm::mat4 camToWorld = cam.invVMatrix();

	//cimg_library::CImg<unsigned char> image(resolution.x, resolution.y, 1, 3, 0);
	for (int y = 0; y < resolution.y; y++) {
		int yind = (resolution.y - 1 - y);
		for (int x = 0; x < resolution.x; x++) {
			const glm::vec2 pixpos { (x + .5f) / resolution.x * 2.0f - 1.0f, 1.0f - (y + .5f) / resolution.y * 2.0f	};
			Ray ray = cam.pixRayDirection(pixpos);//raydirection(resolution, cam, glm::ivec2(x, y), camToWorld, perspective);
			float tri = pixels[3*( yind * resolution.x + x)];
			int tri_id = 0;
			if (tri > 0) {
				// always add ray to tree to calculate size of a node
				tri_id = int(triNum * tri) - 1;
				rst->putPrimitive(ray, tri_id);
				//const unsigned char color[] = { 255 , 0, 0 };
				//image.draw_point(x, y, color);
			}

		}
	}
	//image.save("file.bmp");
}

void cameraSamples(Camera& cam, Shader& rstshader, GLuint rsttex, int no_triangles, 
					glm::ivec2 resolution, std::vector<Vertex>& vertices, RaySpaceTree* rst) {
	glm::mat4 mvp = cam.vpMatrix();
	glUniformMatrix4fv(glGetUniformLocation(rstshader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
	glUniform3fv(glGetUniformLocation(rstshader.index, "viewPos"), 1, glm::value_ptr(cam.position));

	glClearDepth(1.0f);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Execute draw command
	glDrawArrays(GL_TRIANGLES, 0, vertices.size());

	int size = resolution.x * resolution.y;
	glBindTexture(GL_TEXTURE_2D, rsttex);
	GLfloat* pixels = new GLfloat[size * 3];
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pixels);
	makeSample(resolution, cam, pixels, no_triangles, rst);
}

void setupTextureRender(Shader& rstshader, int rstwidth, int rstheight, GLuint vao, GLuint &rsttex, GLuint &framebuffer) {
	rsttex = Buffers::generateTexture(rstwidth, rstheight);
	framebuffer = Buffers::generateFrameBuffer(rsttex);
	GLuint depthbuffer = Buffers::generateDepthBuffer(rstwidth, rstheight);
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) return;

	// Bind the shader
	glUseProgram(rstshader.index);

	// Bind vertex data & buffer
	glBindVertexArray(vao);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

	// Set viewport size
	glViewport(0, 0, rstwidth, rstheight);
	glEnable(GL_DEPTH_TEST);
}

void makeRST(RaySpaceTree* rst, std::vector<Vertex> &vertices, GLuint vao, int no_triangles, 
			Shader &shader, Shader &rstshader, int rstwidth, int rstheight, int option) {

	GLuint rsttex, framebuffer;
	setupTextureRender(rstshader, rstwidth, rstheight, vao, rsttex, framebuffer);

	// Initialize tree
	rst->depth = 8;
	rst->mainDir = glm::vec3(1, 0, 0);
	rst->construct(0, vertices, option);

	// Generate cameras
	// w,h,near,far all depend on bounding sphere!
	// Check if this works!!!
	Orthocamera cam(2, 2, 0.1f, 30.f);
	//PersCamera cam(0.1f, 30.f);

	//need cube to check intersection - if no intersection, not in rst!!
	// need center of bounding sphere and radius of bounding sphere

	///////////////// Sample Sphere implementation instead //////////////////////////
	int noSample = 10;
	float sampleWidth = 2.f / (float)noSample;
	float sampleStart[2] = { 0.f, -1.f };
	for (int z = 0; z < noSample; z++) {
		for (int y = 0; y < noSample; y++) {
			cam.setPosition(glm::vec3(-2.f, sampleStart[0] + sampleWidth*y, sampleStart[0] + sampleWidth * z));
			cam.forward = -cam.position; //should be center of bounding sphere (check if 0,0,0 works)
			cameraSamples(cam, rstshader, rsttex, no_triangles, glm::ivec2(rstwidth, rstheight), vertices, rst);
		}
	}

	//GLuint quadvao = Models::quadGeneration();
	//while (!glfwWindowShouldClose(window))
	//{
	//	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	//	glDisable(GL_DEPTH_TEST); // disable depth test so screen-space quad isn't discarded due to depth test.
	//	// clear all relevant buffers
	//	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // set clear color to white (not really necessary actually, since we won't be able to see behind the quad anyways)
	//	glClear(GL_COLOR_BUFFER_BIT);
	//	glUseProgram(shader);
	//	glBindVertexArray(quadvao);
	//	glBindTexture(GL_TEXTURE_2D, rsttex);	// use the color attachment texture as the texture of the quad plane	
	//	glDrawArrays(GL_TRIANGLES, 0, 6);
	//	glfwSwapBuffers(window);
	//}
	glDeleteFramebuffers(1, &framebuffer);
	glDeleteTextures(1, &rsttex);
}

void getRstStatistics(RaySpaceTree *rst, int no_triangles) {
	// Some statistics
	std::vector<int> nodenr = std::vector<int>(rst->nodes.size());
	std::vector<int> duplicates = rst->countDuplicates(no_triangles, nodenr);
	for (int i = 0; i < duplicates.size(); i++) {
		if (duplicates[i] > 0)
			std::cout << i << " " << duplicates[i] << std::endl;
	}
}

int main() {
	if (!glfwInit()) {
		std::cerr << "Failed to initialize GLFW!" << std::endl;
		return EXIT_FAILURE;
	}

	//////////////////// Create window and OpenGL 4.3 debug context
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


	GLFWwindow* window = glfwCreateWindow(width, height, "Shadow mapping practical", nullptr, nullptr);
	if (!window) {
		std::cerr << "Failed to create OpenGL context!" << std::endl;
		std::cout << "Press enter to close."; getchar();
		return EXIT_FAILURE;
	}

	glfwSetKeyCallback(window, keyboardHandler);
	glfwSetMouseButtonCallback(window, mouseButtonHandler);
	glfwSetCursorPosCallback(window, cursorPosHandler);
	glfwSetScrollCallback(window, scrollHandler);

	// Activate the OpenGL context
	glfwMakeContextCurrent(window);

	// Initialize GLEW extension loader
	glewExperimental = GL_TRUE;
	glewInit();

	// Set up OpenGL debug callback
	glDebugMessageCallback(debugCallback, nullptr);

	////////////////////////// Create shaders
	Shader mainProgram = Shader("standard.vert", "standard.frag");
	Shader rstProgram = Shader("rst.vert", "rst.frag");
	Shader texprojProgram = Shader("texproj.vert", "texproj.frag");
	Shader lineProgram = Shader("line.vert", "line.frag");

	////////////////////////// Load vertices of model
	std::vector<Vertex> vertices;
	int noPrimitives = 0;
	const char *filename = "scene.obj";
	GLuint vao = Models::loadModel(vertices, noPrimitives, filename);

	/////////////////////// Create cube and lines for cube
	float cubeSize = 2.f;
	glm::vec3 cubeMin = glm::vec3(-1.f, 0.f, -1.f);
	glm::vec3 cubeMax = cubeMin + cubeSize;
	Cube cube(cubeMin, cubeSize);
	GLuint lineVAO = Models::cubelineGeneration(cubeMin, cubeMax);

	/////////////////////// Create sphere
	int spherenum = 0;
	GLuint spherevbo;
	GLuint sphereVAO = Models::sphere(10, 20, spherenum, 1.f, spherevbo);

	/////////////////// Create main camera
	PreviewCamera mainCamera;
	mainCamera.aspect = width / (float)height;

	////////////////// Unbind the off-screen framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	/////////////////// Make rayspacetree
	std::vector<glm::vec3> splitters;
	std::vector<glm::vec3> leafRays;
	std::vector<glm::vec3> rayColors;

	const int w = 400, h = 400;
	RaySpaceTree rst(&cube);
	makeRST(&rst, vertices, vao, noPrimitives, texprojProgram, rstProgram, w, h, 0);
	// Use rst info to visualize
	splitters = rst.getSplittingLinesInCube();

	GLuint splitterVao = Models::vaoLineGeneration(splitters);
	GLuint leafRayVao;
	int noLeaves = std::pow(2, rst.depth);
	int leafnum = leafNumber;
	bool showing = false;

	//getRstStatistics(&rst, noPrimitives);

	// Main loop
	glEnable(GL_DEPTH_TEST);

	while (!glfwWindowShouldClose(window)) {

		glfwPollEvents();

		if (leafnum != leafNumber) {
			leafRays = std::vector<glm::vec3>();
			splitters = std::vector<glm::vec3>();
			rayColors = std::vector<glm::vec3>();

			while (leafRays.size() == 0) {
				leafNumber = leafNumber % noLeaves;
				rst.getViewingLinesInLeaf(leafNumber, leafRays, rayColors, splitters);
				splitterVao = Models::vaoLineGeneration(splitters);
				leafNumber++;
			}

			leafRayVao = Models::vaoLineGenerationWithColor(leafRays, rayColors, &cube);
			showing = true;

			leafnum = leafNumber;
		}

		updateCamera(mainCamera, width, height);
		//glm::mat4 mvp = mainCamera.vpmMatrix();
		glm::mat4 mvp = mainCamera.vpmMatrix();

		// Bind the shader
		glUseProgram(mainProgram.index);
		glUniformMatrix4fv(glGetUniformLocation(mainProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

		// Set viewport size
		glViewport(0, 0, width, height);

		// Clear the framebuffer to black and depth to maximum value
		glClearDepth(1.0f);  
		glClearColor(0.1f, 0.2f, 0.3f, 1.0f); 
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glBindVertexArray(vao);
		glDrawArrays(GL_TRIANGLES, 0, vertices.size());

		glBindVertexArray(lineVAO);
		glDrawArrays(GL_LINES, 0, 24);

		glBindVertexArray(splitterVao);
		glDrawArrays(GL_LINES, 0, splitters.size());

		//glBindVertexArray(sphereVAO);
		//glEnable(GL_PRIMITIVE_RESTART);
		//glPrimitiveRestartIndex(GL_PRIMITIVE_RESTART_FIXED_INDEX);
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spherevbo);
		//glDrawElements(GL_QUAD_STRIP, spherenum, GL_UNSIGNED_INT, NULL);

		if (showing) {
			glUseProgram(lineProgram.index);
			glUniformMatrix4fv(glGetUniformLocation(lineProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			glBindVertexArray(leafRayVao);
			glDrawArrays(GL_LINES, 0, leafRays.size());
		}

		// Present result to the screen
		glfwSwapBuffers(window);
	}
	
	glfwDestroyWindow(window);
	
	glfwTerminate();

    return 0;
}