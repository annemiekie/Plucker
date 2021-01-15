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
#include "raytracer.h"
#include "model.h"
#include "sphere.h"
#include "sphereSampler.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>

#include <CImg.h>


// Configuration
const int width = 800;
const int height = 800;
bool update = false;
bool trace = true;
bool togglePoints = false;
bool toggleSphere = false;
bool toggleBbox = false;
bool toggleLines = false;
bool drawRay = false;
glm::vec2 drawRayPos(0, 0);

// Key handle function
void keyboardHandler(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	//cameraKeyboardHandler(key, action);
	if (action == GLFW_PRESS) {

		switch (key)
		{
		case GLFW_KEY_1:
			update = true;
			break;
		case GLFW_KEY_2:
			togglePoints = !togglePoints;
			break;
		case GLFW_KEY_3:
			toggleSphere = !toggleSphere;
			break;
		case GLFW_KEY_4:
			toggleBbox = !toggleBbox;
			break;
		case GLFW_KEY_5:
			toggleLines = !toggleLines;
			break;
		default:
			break;
		}
	}
}

// Mouse button handle function
void mouseButtonHandler(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_2 && action == GLFW_RELEASE) {
		double xpos, ypos;
		//getting cursor position
		glfwGetCursorPos(window, &xpos, &ypos);
		drawRay = true;
		drawRayPos = { (xpos + .5f) / float(width) * 2.0f - 1.0f, 1.0f - (ypos + .5f) / float(height) * 2.0f };

	}
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

void makeSample(glm::ivec2 resolution, Camera& cam, GLfloat* pixels, RaySpaceTree* rst) {

	for (int y = 0; y < resolution.y; y++) {
		int yind = (resolution.y - 1 - y);
		for (int x = 0; x < resolution.x; x++) {
			const glm::vec2 pixpos { (x + .5f) / resolution.x * 2.0f - 1.0f, 1.0f - (y + .5f) / resolution.y * 2.0f	};
			Ray ray = cam.pixRayDirection(pixpos);
			float tri = pixels[yind * resolution.x + x];
			if (tri > 0) {
				// always add ray to tree to calculate size of a node
				int tri_id = int(tri - 1);
				rst->putPrimitive(ray, tri_id);
			}
		}
	}
}

void makeSample(glm::ivec2 res, Camera& cam, GLfloat* pixels, std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int& count, int& zerocount) {

	for (int y = 0; y < res.y; y++) {
		int yind = (res.y - 1 - y);
		for (int x = 0; x < res.x; x++) {
			const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
			Ray ray = cam.pixRayDirection(pixpos);
			float tri = pixels[yind * res.x + x];
			samples.push_back({ count, int(tri - 1) });
			count++;
			if (tri > 0) tris.insert(int(tri - 1));
			else zerocount++;
		}
	}
}

GLfloat* cameraSamples(Camera& cam, Shader& rstshader, GLuint rsttex, Model& model, glm::ivec2 resolution) {
	glm::mat4 mvp = cam.vpMatrix();
	glUniformMatrix4fv(glGetUniformLocation(rstshader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
	glUniform3fv(glGetUniformLocation(rstshader.index, "viewPos"), 1, glm::value_ptr(cam.position));

	glClearDepth(1.0f);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Execute draw command
	glDrawArrays(GL_TRIANGLES, 0, model.vertices.size());

	int size = resolution.x * resolution.y;
	glBindTexture(GL_TEXTURE_2D, rsttex);
	GLfloat* pixels = new GLfloat[size];
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, pixels);
	return pixels;
}

void setupTextureRender(Shader& rstshader, int rstwidth, int rstheight, GLuint vao, GLuint &rsttex, GLuint &framebuffer) {

	glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);

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
	glEnable(GL_CULL_FACE);
}

void makeRST(RaySpaceTree* rst, Model& model, SphereSampler& sampler, 
			 Shader &rstshader, int rstwidth, int rstheight, int option) {

	GLuint rsttex, framebuffer;

	setupTextureRender(rstshader, rstwidth, rstheight, model.vao, rsttex, framebuffer);

	// Initialize tree
	rst->construct(model.vertices, option);

	// Generate cameras
	Orthocamera cam(model.radius, model.radius, 0.1f, 30.f);
	glm::ivec2 resolution(rstwidth, rstheight);

	//need cube to check intersection - if no intersection, not in rst!!
	for (const glm::vec3 &sample : sampler.samples) {
		cam.setPositionAndForward(sample, model.center);
		GLfloat *pixels = cameraSamples(cam, rstshader, rsttex, model, resolution);
		makeSample(resolution, cam, pixels, rst);
	}
	glDeleteFramebuffers(1, &framebuffer);
	glDeleteTextures(1, &rsttex);
}

void makeAdaptiveRST(RaySpaceTree* rst, Model& model, SphereSampler& sampler,
					 Shader& rstshader, int rstwidth, int rstheight, int option) {

	GLuint rsttex, framebuffer;

	setupTextureRender(rstshader, rstwidth, rstheight, model.vao, rsttex, framebuffer);

	// Generate cameras
	Orthocamera cam(model.radius, model.radius, 0.1f, 30.f);
	glm::ivec2 resolution(rstwidth, rstheight);

	int imgsize = rstwidth * rstheight;
	std::vector<std::pair<int, int>> samples = std::vector<std::pair<int,int>>(); // or sample(int raynr, int trinr)
	std::set<int> tris = std::set<int>();
	std::vector<Orthocamera> cams = std::vector<Orthocamera>();
	int count = 0;
	int zerocount = 0;

	std::cout << "Rasterizing camera samples..." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	// prepare samples
	for (int i = 0; i < sampler.samples.size(); i++) {// const glm::vec3 & sample : sampler.samples) {
		cam.setPositionAndForward(sampler.samples[i], model.center);
		GLfloat* pixels = cameraSamples(cam, rstshader, rsttex, model, resolution);
		makeSample(resolution, cam, pixels, samples, tris, count, zerocount);
		cams.push_back(cam);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Made and stored " << rstwidth*rstheight*sampler.samples.size() << " samples in " << diff << " ms." << std::endl;

	std::cout << "Constructing the tree..." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();
	rst->constructAdaptive(model.vertices, resolution, cams, samples, tris, samples.size(), zerocount);
	end_time = std::chrono::high_resolution_clock::now();
	diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Constructed RST in " << diff << " ms." << std::endl;

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
	std::cout << std::endl;
	for (int i = 0; i < nodenr.size(); i++) {
		if (nodenr[i] > 0)
			std::cout << i << " " << nodenr[i] << std::endl;
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
	std::cout << "Set up shader programs." << std::endl;


	////////////////////////// Load vertices of model
	const char *filename = "scene.obj";
	Model model(filename, 1);

	/////////////////////// Create cube and lines for cube
	float cubeSize = 2.f;
	//glm::vec3 cubeMin = glm::vec3(-1.f, 0.f, -1.f);
	//glm::vec3 cubeMax = cubeMin + cubeSize;
	Cube cube(model.center, cubeSize);
	GLuint cubeVAO = cube.cubelineGeneration();
	std::cout << "Generated viewing cube." << std::endl;

	////////////////////// Create show bounding box
	GLuint bboxVAO = model.boundingBox.cubelineGeneration();

	/////////////////////// Create show sphere
	Sphere sphere(model.center, model.radius);
	sphere.vaoGeneration(10, 20);

	///////////////////// Create sample locations
	Sphere sampleSphere(model.center, 2.f);
	sampleSphere.vaoGeneration(10, 20);
	SphereSampler sampler(&sampleSphere);
	int noSamples = 100;
	char dir = 'X';
	sampler.createSamples(-1, dir, noSamples);
	sampler.vaoGeneration();
	glm::vec3 mainDir(1, 0, 0); /// check thisssss!!!!!!!!!!!
	std::cout << "Created " << sampler.samples.size() << " camera locations in " << dir << " direction" << std::endl;

	/////////////////// Create main camera
	PreviewCamera mainCamera;
	mainCamera.aspect = width / (float)height;

	////////////////// Unbind the off-screen framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	/////////////////// Make rayspacetree
	const int w = 400, h = 400;
	RaySpaceTree rst(&cube);
	rst.depth = 8;
	rst.mainDir = glm::vec3(1, 0, 0);
	std::cout << "Constructing RST with depth = " << rst.depth << " and camera sample size = " << w << "x" << h << "..." << std::endl;
	makeAdaptiveRST(&rst, model, sampler, rstProgram, w, h, 0);
	//makeRST(&rst, model, sampler, rstProgram, w, h, 0);
	rst.printTree();

	////////////////// Check result by RayTracing
	if (trace) {
		std::cout << "Raytracing to check tree..." << std::endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		Orthocamera raytracecam(model.radius, model.radius, 0.1f, 30.f);
		raytracecam.setPositionAndForward(glm::vec3(-2.f, 1.f, 0.f), model.center);
		RayTracer::imageTracer(raytracecam, &rst, glm::ivec2(100, 100), model);
		auto end_time = std::chrono::high_resolution_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
		std::cout << "Raytraced image in " << diff << " ms." << std::endl;
	}

	////////////////// Use rst info to visualize
	std::vector<glm::vec3> splitters = rst.getSplittingLinesInCube();;
	std::vector<glm::vec3> leafRays;
	std::vector<glm::vec3> rayColors;
	GLuint splitterVao = Models::vaoLineGeneration(splitters);
	GLuint leafRayVao = 0;
	int noLeaves = std::pow(2, rst.depth);
	int leafnum = -1;
	bool showing = false;
	bool viewline = false;
	float alphaLines = 1.f;
	GLuint rayVAO;

	//getRstStatistics(&rst, model.primsize);

	// Main loop
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_CCW);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	while (!glfwWindowShouldClose(window)) {

		glfwPollEvents();
		
		if (update) {
			alphaLines = 1.f;
			leafRays = std::vector<glm::vec3>();
			rayColors = std::vector<glm::vec3>();

			while (leafRays.size() == 0) {
				splitters = std::vector<glm::vec3>();
				leafnum++;
				leafnum = leafnum % noLeaves;
				rst.getViewingLinesInLeaf(leafnum, leafRays, rayColors, splitters);
			}
			splitterVao = Models::vaoLineGeneration(splitters);
			leafRayVao = Models::vaoLineGenerationWithColor(leafRays, rayColors, &cube);
			toggleLines = true;
			update = false;
			viewline = false;
		}

		if (drawRay) {
			leafRays = std::vector<glm::vec3>();
			splitters = std::vector<glm::vec3>();
			rayColors = std::vector<glm::vec3>();

			alphaLines = 0.5f;
			Ray r = mainCamera.pixRayDirection(drawRayPos);
			rst.descendWithLines(r, splitters, leafRays, rayColors);
			glm::vec3 color(0,1,0);
			rayVAO = Models::rayVao(r, color, &cube, mainDir);

			splitterVao = Models::vaoLineGeneration(splitters);
			leafRayVao = Models::vaoLineGenerationWithColor(leafRays, rayColors, &cube);

			drawRay = false;
			viewline = true;
		}

		updateCamera(mainCamera, width, height);
		glm::mat4 mvp = mainCamera.vpMatrix();

		// Bind the shader
		glUseProgram(mainProgram.index);
		glUniformMatrix4fv(glGetUniformLocation(mainProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

		// Set viewport size
		glViewport(0, 0, width, height);

		// Clear the framebuffer to black and depth to maximum value
		glClearDepth(1.0f);  
		glClearColor(0.1f, 0.2f, 0.3f, 1.0f); 
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glBindVertexArray(model.vao);
		glDrawArrays(GL_TRIANGLES, 0, model.vertices.size());

		glBindVertexArray(cubeVAO);
		glDrawArrays(GL_LINES, 0, 24);

		if (toggleBbox) {
			glBindVertexArray(bboxVAO);
			glDrawArrays(GL_LINES, 0, 24);
		}

		if (togglePoints) {
			glBindVertexArray(sampler.vao);
			glDrawArrays(GL_POINTS, 0, sampler.samples.size());
		}

		glUseProgram(lineProgram.index);
		glUniformMatrix4fv(glGetUniformLocation(lineProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
		glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), (GLfloat)alphaLines);

		if (toggleLines) {
			glBindVertexArray(splitterVao);
			glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), 1.f);
			glDrawArrays(GL_LINES, 0, splitters.size());
			glBindVertexArray(leafRayVao);
			glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), (GLfloat)alphaLines);
			glDrawArrays(GL_LINES, 0, leafRays.size());
		}

		if (viewline) {
			glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), 1.f);
			glBindVertexArray(rayVAO);
			glDrawArrays(GL_LINES, 0, 2);
		}

		if (toggleSphere) {
			glUseProgram(mainProgram.index);
			mvp = mainCamera.vpmMatrix(sphere.center);
			glUniformMatrix4fv(glGetUniformLocation(mainProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			glBindVertexArray(sphere.vao);
			glDrawElements(GL_LINES, sphere.indSize, GL_UNSIGNED_INT, (void*)0);

			glBindVertexArray(sampleSphere.vao);
			glDrawElements(GL_LINES, sampleSphere.indSize, GL_UNSIGNED_INT, (void*)0);
		}

		// Present result to the screen
		glfwSwapBuffers(window);
	}
	
	glfwDestroyWindow(window);
	
	glfwTerminate();

    return 0;
}