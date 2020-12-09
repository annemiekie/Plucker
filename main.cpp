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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <CImg.h>


// Configuration
const int width = 800;
const int height = 800;
bool update = false;
bool trace = false;

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

	//cimg_library::CImg<unsigned char> image(resolution.x, resolution.y, 1, 3, 0);
	for (int y = 0; y < resolution.y; y++) {
		int yind = (resolution.y - 1 - y);
		for (int x = 0; x < resolution.x; x++) {
			const glm::vec2 pixpos { (x + .5f) / resolution.x * 2.0f - 1.0f, 1.0f - (y + .5f) / resolution.y * 2.0f	};
			Ray ray = cam.pixRayDirection(pixpos);
			float tri = pixels[yind * resolution.x + x];
			int tri_id = 0;
			if (tri > 0) {
				// always add ray to tree to calculate size of a node
				tri_id = int(tri - 1);//int(triNum * tri - 1.f);
				rst->putPrimitive(ray, tri_id);
				//const unsigned char color[] = { 255 , 0, 0 };
				//image.draw_point(x, y, color);
			}

		}
	}
	//image.save("file.bmp");
}

void cameraSamples(Camera& cam, Shader& rstshader, GLuint rsttex, Model& model, 
					glm::ivec2 resolution, RaySpaceTree* rst) {
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
	makeSample(resolution, cam, pixels, model.primsize, rst);
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

void makeRST(RaySpaceTree* rst, Model& model, Shader &shader, Shader &rstshader, int rstwidth, int rstheight, int option, int noSample) {

	GLuint rsttex, framebuffer;
	setupTextureRender(rstshader, rstwidth, rstheight, model.vao, rsttex, framebuffer);

	// Initialize tree
	rst->depth = 8;
	rst->mainDir = glm::vec3(1, 0, 0);
	rst->construct(0, model.vertices, option);

	// Generate cameras
	// w,h,near,far all depend on bounding sphere!
	Orthocamera cam(model.radius, model.radius, 0.1f, 30.f);
	//PersCamera cam(0.1f, 30.f);

	//need cube to check intersection - if no intersection, not in rst!!

	///////////////// Sample Sphere implementation instead //////////////////////////
	float sampleWidth = 2.f / (float)(noSample - 1);
	float sampleStart[2] = { 0.001f, -1.001f };

	//#pragma omp parallel for
	for (int z = 0; z < noSample; z++) {
		for (int y = 0; y < noSample; y++) {
			//cam.setPositionAndForward(glm::vec3(-2.f, sampleStart[0] + sampleWidth*y, sampleStart[1] + sampleWidth * z));
			glm::vec3 position = glm::vec3(-2.f, sampleStart[0] + sampleWidth * y, sampleStart[1] + sampleWidth * z); //replace with sample
			cam.setPositionAndForward(position, model.center);
			cameraSamples(cam, rstshader, rsttex, model, glm::ivec2(rstwidth, rstheight), rst);
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

GLuint makePoints(int &size) {
	std::vector<GLfloat> points;
	float phi, theta, sintheta;
	float x, y, z;
	float factorPhi = 2.f * glm::pi<float>() * 2.f / (1.f + sqrtf(5));
	int N = 1000;
	float factorTheta = 2.f / (2.f * N + 1.f);
	float r = 1.5f;

	float maxx = -1000.f, maxy = -1000.f, maxz = -1000.f;

	for (int i = -N; i <= N; i++) {
		phi = i * factorPhi;
		sintheta = i * factorTheta;
		theta = asin(sintheta);
		x = r * cos(theta) * cos(phi);
		y = r * cos(theta) * sin(phi);
		z = r * sintheta;

		points.push_back(x);
		points.push_back(y);
		points.push_back(z);
	}

	GLuint pointVAO, pointVBO;
	glGenVertexArrays(1, &pointVAO);
	glBindVertexArray(pointVAO);

	glGenBuffers(1, &pointVBO);
	glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*points.size(), &points[0], GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*3, (void*)0);
	size = points.size()/3;
	return pointVAO;

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
	const char *filename = "scene.obj";
	Model model(filename, 1);

	/////////////////////// Create cube and lines for cube
	float cubeSize = 2.f;
	glm::vec3 cubeMin = glm::vec3(-1.f, 0.f, -1.f);
	glm::vec3 cubeMax = cubeMin + cubeSize;
	Cube cube(cubeMin, cubeSize);
	GLuint cubeVAO = cube.cubelineGeneration();

	////////////////////// Create show bounding box
	GLuint bboxVAO = model.boundingBox.cubelineGeneration();

	/////////////////////// Create show sphere
	Sphere sphere(model.center, model.radius);
	sphere.vaoGeneration(10, 20);

	///////////////////// Create sample locations show
	int pointsize = 0;
	GLuint pointsVAO = makePoints(pointsize);

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
	int noSample = 10;
	makeRST(&rst, model, texprojProgram, rstProgram, w, h, 0, noSample);

	////////////////// Check result by RayTracing
	Orthocamera raytracecam(model.radius, model.radius, 0.1f, 30.f);
	raytracecam.setPositionAndForward(glm::vec3(-2.f, 1.f, 0.f), model.center);
	if (trace) RayTracer::imageTracer(raytracecam, &rst, glm::ivec2(w,h), model);

	////////////////// Use rst info to visualize
	splitters = rst.getSplittingLinesInCube();
	GLuint splitterVao = Models::vaoLineGeneration(splitters);
	GLuint leafRayVao;
	int noLeaves = std::pow(2, rst.depth);
	int leafnum = -1;
	bool showing = false;

	//getRstStatistics(&rst, noPrimitives);

	// Main loop
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_CCW);
	glEnable(GL_PROGRAM_POINT_SIZE);

	while (!glfwWindowShouldClose(window)) {

		glfwPollEvents();
		
		if (update) {
			leafRays = std::vector<glm::vec3>();
			splitters = std::vector<glm::vec3>();
			rayColors = std::vector<glm::vec3>();

			while (leafRays.size() == 0) {
				leafnum++;
				leafnum = leafnum % noLeaves;
				rst.getViewingLinesInLeaf(leafnum, leafRays, rayColors, splitters);
				splitterVao = Models::vaoLineGeneration(splitters);
			}

			leafRayVao = Models::vaoLineGenerationWithColor(leafRays, rayColors, &cube);
			showing = true;
			update = false;
		}

		updateCamera(mainCamera, width, height);
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

		glBindVertexArray(model.vao);
		glDrawArrays(GL_TRIANGLES, 0, model.vertices.size());

		glBindVertexArray(cubeVAO);
		glDrawArrays(GL_LINES, 0, 24);

		glBindVertexArray(bboxVAO);
		glDrawArrays(GL_LINES, 0, 24);

		glBindVertexArray(splitterVao);
		glDrawArrays(GL_LINES, 0, splitters.size());

		if (showing) {
			glUseProgram(lineProgram.index);
			glUniformMatrix4fv(glGetUniformLocation(lineProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			glBindVertexArray(leafRayVao);
			glDrawArrays(GL_LINES, 0, leafRays.size());
		}
		glUseProgram(mainProgram.index);
		mvp = mainCamera.vpmMatrix(sphere.center);
		glUniformMatrix4fv(glGetUniformLocation(mainProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
		glBindVertexArray(sphere.vao);
		glDrawElements(GL_LINES, sphere.indSize, GL_UNSIGNED_INT, (void*)0);

		glBindVertexArray(pointsVAO);
		glDrawArrays(GL_POINTS, 0, pointsize);

		// Present result to the screen
		glfwSwapBuffers(window);
	}
	
	glfwDestroyWindow(window);
	
	glfwTerminate();

    return 0;
}