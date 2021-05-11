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
#include "geoObject.h"
#include "linemodel.h"

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
bool trace = false;
bool togglePoints = false;
bool toggle4lines = false;
bool toggleSphere = false;
bool toggleBbox = false;
bool toggleLines = true;
bool toggleSampleLines = false;
bool nextleaf = false;
bool drawRay = false;
bool print = false;
glm::vec2 drawRayPos(0, 0);
bool togglePicking = false;
glm::vec2 pickingPos(0, 0);
bool picking = false;
bool sphereOrCube = false;

// OpenGL debug callback
void APIENTRY debugCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {
	if (severity != GL_DEBUG_SEVERITY_NOTIFICATION) {
		std::cerr << "OpenGL: " << message << std::endl;
	}
}


bool checkTriIntersect(int i, Ray& r, RaySpaceTree* rst) {

	Ray edge;
	edge = Ray(rst->model->vertices[i + 1].pos, rst->model->vertices[i].pos);
	if (edge.side(r)) return false;
	edge = Ray(rst->model->vertices[i + 2].pos, rst->model->vertices[i + 1].pos);
	if (edge.side(r)) return false;
	edge = Ray(rst->model->vertices[i].pos, rst->model->vertices[i + 2].pos);
	if (edge.side(r)) return false;
	return true;
}

void makeSample(glm::ivec2 res, Camera* cam, GLfloat* pixels, RaySpaceTree* rst, bool storeRays, 
				std::vector<Ray>& raytraceRays) {

	for (int y = 0; y < res.y; y++) {
		int yind = (res.y - 1 - y);
		for (int x = 0; x < res.x; x++) {
			const glm::vec2 pixpos { (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
			Ray ray = cam->pixRayDirection(pixpos);
			float tri = pixels[yind * res.x + x];
			if (tri > 0) {
				int tri_id = int(tri - 1);
				if (checkTriIntersect(tri_id*3, ray, rst)) 
					rst->putPrimitive(ray, tri_id, storeRays);
				else raytraceRays.push_back(ray);

			}
		}
	}
}

int traceMissingRay(Ray ray, RaySpaceTree *rst) {
	float depth = 1000;
	int prim = -1;
	float diff;
	for (int tri_id = 0; tri_id < rst->model->vertices.size() / 3; tri_id++) {
		float t = RayTracer::intersection(tri_id * 3, rst->model->vertices, ray, diff);
		if (t > 0 && t < depth) {
			prim = tri_id;
		}
	}
	return prim;
}

void makeSample(glm::ivec2 res, Camera* cam, GLfloat* pixels, RaySpaceTree* rst, std::vector<std::pair<int, int>>& samples,
				std::set<int>& tris, int& count) {

	for (int y = 0; y < res.y; y++) {
		int yind = (res.y - 1 - y);
		for (int x = 0; x < res.x; x++) {
			const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
			Ray ray = cam->pixRayDirection(pixpos);
			float tri = pixels[yind * res.x + x];
			int tri_id = int(tri - 1);
			if (tri_id > -2) {
				if (tri_id >= 0) {
					if (checkTriIntersect(tri_id * 3, ray, rst)) tris.insert(tri_id);
					else {
						tri_id = traceMissingRay(ray, rst);
						if (tri_id >= 0) tris.insert(tri_id);
					}
					samples.push_back({ count, tri_id});
				}
				else if (tri_id == -1) samples.push_back({ count, -1 });
			}
			count++;
		}
	}
}

GLfloat* cameraSamples(Camera* cam, Shader& rstshader, GLuint rsttex, Model* model, glm::ivec2 resolution, int second) {
	glm::mat4 mvp = cam->vpMatrix();
	glUniformMatrix4fv(glGetUniformLocation(rstshader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
	//glUniform1i(glGetUniformLocation(rstshader.index, "secondpass"), second);
		//GLuint texUnitLoc = glGetUniformLocation(rsttex, "texUnit");
		//glProgramUniform1i(rsttex, texUnitLoc, 0);

	glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, secondtex);
	glBindTexture(GL_TEXTURE_2D, rsttex);

	glClearDepth(1.0f);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Execute draw command
	glDrawArrays(GL_TRIANGLES, 0, model->vertices.size());


	int size = resolution.x * resolution.y;
	GLfloat* pixels = new GLfloat[size];
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, pixels);

	return pixels;
}

void setupTextureRender(Shader& rstshader, int rstwidth, int rstheight, GLuint vao, GLuint &rsttex, GLuint &depthtex, GLuint &framebuffer) {

	glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
	rsttex = Buffers::generateTexture(rstwidth, rstheight);
	depthtex = Buffers::generateDepthTexture(rstwidth, rstheight);
	framebuffer = Buffers::generateFrameBuffer(rsttex, depthtex);
	//GLuint depthbuffer = Buffers::generateDepthBuffer(rstwidth, rstheight);
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) return;


	// Bind the shader
	glUseProgram(rstshader.index);

	// Bind vertex data & buffer
	glBindVertexArray(vao);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

	// Set viewport size
	glViewport(0, 0, rstwidth, rstheight);
	glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LESS);
	glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
}

void fillMoreSamples(RaySpaceTree* rst, SphereSampler& sampler,
					Shader& rstshader, GLuint rsttex, glm::ivec2 res, Camera* cam, bool storeRays) {

	std::cout << "Rasterizing camera samples to fill tree..." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	std::vector<Ray> raytraceRays;
	for (int i = 0; i < sampler.samples.size(); i++) {
		if (i % sampler.ratio == 0) continue;
		cam->setPositionAndForward(sampler.samples[i], rst->model->center);
		GLfloat* pixels = cameraSamples(cam, rstshader, rsttex, rst->model, res, 0);
		makeSample(res, cam, pixels, rst, false, raytraceRays);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	int noSamples = (sampler.samples.size() - (sampler.samples.size() / sampler.ratio + 1)) * res.x*res.y;
	std::cout << "Put " << noSamples << " samples in RST in " << diff << " ms." << std::endl;

}

Ray makeRayHorizontal(int num, Model* model) {
	Vertex v = model->vertices[num];
	return 	Ray(glm::vec3(v.pos.x, 0, v.pos.z), glm::vec3(v.pos.x, 1,v.pos.z));
}

Ray makeRayVertical(int num, Model* model) {
	Vertex v = model->vertices[num];
	return 	Ray(glm::vec3(v.pos.x, v.pos.y, 0), glm::vec3(v.pos.x, v.pos.y, 1));
}

std::vector<Ray> constructRays(Model* model) {
	std::vector<Ray> rays = std::vector<Ray>();
	rays.push_back(makeRayHorizontal(0, model));
	rays.push_back(makeRayVertical(1000, model));
	rays.push_back(makeRayHorizontal(50000, model));
	rays.push_back(makeRayVertical(400, model));
	rays.push_back(makeRayHorizontal(800, model));
	rays.push_back(makeRayVertical(20000, model));
	rays.push_back(makeRayHorizontal(80, model));
	rays.push_back(makeRayVertical(300, model));
	return rays;
}

std::vector<Ray> constructRaysForPyramid(Model* model) {
	std::vector<Ray> rays = std::vector<Ray>();
	rays.push_back(Ray(model->vertices2[0], model->vertices2[1]));
	rays.push_back(Ray(model->vertices2[2], model->vertices2[1]));
	rays.push_back(Ray(model->vertices2[0], model->vertices2[2]));
	rays.push_back(Ray(model->vertices2[0], model->vertices2[3]));
	rays.push_back(Ray(model->vertices2[1], model->vertices2[3]));
	rays.push_back(Ray(model->vertices2[2], model->vertices2[3]));
	return rays;
}

std::vector<Ray> constructRaysRandom(Model* model, int level) {
	std::vector<Ray> rays = std::vector<Ray>();
	int i = 0;
	while (i < level) {
		float r = (float)rand() / static_cast<float> (RAND_MAX);
		int r1 = int(r * (model->vertices.size() / 3.f));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);

		Ray ray(model->vertices[3 * r1 + redge].pos, model->vertices[3 * r1 + (redge+1)%3].pos);
		bool dupli = false;
		for (auto r : rays) if (r.equal(ray, 1E-6)) { dupli = true; break; }
		if (dupli) continue;
		else {
			rays.push_back(ray);
			i++;
		}
	}
	return rays;
}

void makeRST(RaySpaceTree* rst, SphereSampler& sampler, Shader &rstshader, glm::ivec2 res, 
			Camera* cam, GLuint rsttex, int option, bool storeRays) {

	std::vector<Ray> rays = constructRaysRandom(rst->model, rst->depth);
	// Initialize tree
	rst->construct(option, rays);

	// Generate cameras
	std::vector<Ray> raytraceRays;

	std::cout << "Rasterizing camera samples to create tree..." << std::endl;

	for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
		std::cout << i << std::endl;
		cam->setPositionAndForward(sampler.samples[i], rst->model->center);
		GLfloat *pixels = cameraSamples(cam, rstshader, rsttex, rst->model, res, 0);
		makeSample(res, cam, pixels, rst, storeRays, raytraceRays);
	}

	//std::cout << "Raytracing " << raytraceRays.size() << " faulty rays..." << std::endl;

	//for (Ray ray : raytraceRays) {
	//	float depth = 1000;
	//	int prim = -1;
	//	float diff;
	//	for (int tri_id = 0; tri_id < rst->model->vertices.size() / 3; tri_id++) {
	//		float t = RayTracer::intersection(tri_id * 3, rst->model->vertices, ray, diff);
	//		if (t > 0 && t < depth) {
	//			prim = tri_id;
	//		}
	//	}
	//	if (prim > -1) rst->putPrimitive(ray, prim, storeRays);
	//}
}

void makeAdaptiveRST(RaySpaceTree* rst, SphereSampler& sampler, Shader& rstshader,  
					 glm::ivec2 res, Camera* cam, GLuint rsttex, int option) {

	int imgsize = res.x * res.y;
	std::vector<std::pair<int, int>> samples = std::vector<std::pair<int,int>>();
	std::set<int> tris = std::set<int>();
	std::vector<Camera*> cams = std::vector<Camera*>();
	int count = 0;
	int zerocount = 0;

	std::cout << "Rasterizing camera samples to create tree..." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	std::vector<Ray> raytraceRays;
	// prepare samples
	for (int i = 0; i < sampler.samples.size(); i+=sampler.ratio) {
		Camera* cam1 = cam->makeCopy();
		cam1->setPositionAndForward(sampler.samples[i], rst->model->center);
		
		GLfloat* pixels = cameraSamples(cam1, rstshader, rsttex, rst->model, res, false);
		makeSample(res, cam1, pixels, rst, samples, tris, count);
		cams.push_back(cam);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Made and stored " << samples.size() << " samples in " << diff << " ms." << std::endl;
	std::cout << "Removed " << res.x * res.y * (sampler.samples.size() / sampler.ratio + 1) - samples.size() << " samples for being too close to a primitive edge." << std::endl;

	std::cout << "Constructing the tree..." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();
	if (option == 0) rst->constructAdaptive(res, cams, samples, tris, print);
	else if (option == 1) rst->constructSmartRandom(res, cams, samples, tris);
	end_time = std::chrono::high_resolution_clock::now();
	diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Constructed RST in " << diff << " ms." << std::endl;
}

void checkTreeByTracing(Camera* cam, glm::ivec2 traceRes, RaySpaceTree* rst, Shader rstProgram, GLuint rsttex) {
	std::cout << "Raytracing to check tree..." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();

	std::vector<int> raytraceTri = RayTracer::imageTracer(cam, rst, traceRes);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Raytraced image in " << diff << " ms." << std::endl;

	/////////////// compare with raytrace
	GLfloat* pixels = cameraSamples(cam, rstProgram, rsttex, rst->model, traceRes, false);
	cimg_library::CImg<unsigned char> image(traceRes.x, traceRes.y, 1, 3, 0);
	const unsigned char white[] = { 255 , 255, 255 };
	const unsigned char black[] = { 0 , 0, 0 };

	for (int y = 0; y < traceRes.y; y++) {
		int yind = (traceRes.y - 1 - y);
		for (int x = 0; x < traceRes.x; x++) {
			int triRast = int(pixels[yind * traceRes.x + x]);
			int triTrace = raytraceTri[y * traceRes.x + x];
			if (triRast - triTrace != 0) {
				image.draw_point(x, y, white);
				//std::cout << std::endl << triRast << " " << triTrace << " " << x << " " << y << std::endl;
			}
			//else std::cout << "good ";
			//if (triRast > 0) image.draw_point(x, y, color);

		}
	}
	image.save("rasterizer.bmp");
}

void getRstStatistics(RaySpaceTree* rst, int no_triangles) {
	// Some statistics
	std::vector<int> nodenr = std::vector<int>(rst->nodes.size());
	//std::vector<int> duplicates = rst->countDuplicates(no_triangles, nodenr);
	//std::cout << "DUPLICATES:" << std::endl;
	//for (int i = 0; i < duplicates.size(); i++) {
	//	if (duplicates[i] > 1)
	//		std::cout << i << "'s duplicates: " << duplicates[i] << std::endl;
	//}
	//std::cout << std::endl;
	//for (int i = 0; i < nodenr.size(); i++) {
	//	if (nodenr[i] > 0)
	//		std::cout << i << " " << nodenr[i] << std::endl;
	//}
	int filledleaves = 0;
	int leavesWith0 = 0;
	int leavesWith1 = 0;
	int leavesWithMore = 0;
	int maxdepth = 0;
	int mindepth = 1000;
	for (auto node : rst->nodes) {
		if (node->leaf) {
			filledleaves++;
			if (node->depth > maxdepth) maxdepth = node->depth;
			else if (node->depth < mindepth) mindepth = node->depth;
			if (node->primitiveSet.size() == 0) leavesWith0++;
			else if (node->primitiveSet.size() == 1) leavesWith1++;
			else leavesWithMore++;
		}
	}
	std::cout << "Filled Leaves: " << filledleaves << std::endl;

	std::cout << "Tree Depth Max " << maxdepth << " and Min " << mindepth << std::endl;

	std::cout << "Leaves with 0: " << leavesWith0 << " 1: " << leavesWith1 << " more: " << leavesWithMore << std::endl;


}

// Key handle function
void keyboardHandler(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	//cameraKeyboardHandler(key, action);
	if (action == GLFW_PRESS) {

		switch (key)
		{
		case GLFW_KEY_1:
			nextleaf = true;
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
			toggleSampleLines = !toggleSampleLines;
			update = true;
			break;
		case GLFW_KEY_6:
			toggle4lines = !toggle4lines;
			update = true;
			break;
		case GLFW_KEY_7:
			togglePicking = !togglePicking;
			break;
		case GLFW_KEY_8:
			sphereOrCube = !sphereOrCube;
			update = true;
			break;
		case GLFW_KEY_9:
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
	else if (button == GLFW_MOUSE_BUTTON_1 && action == GLFW_RELEASE) {
		double xpos, ypos;
		//getting cursor position
		glfwGetCursorPos(window, &xpos, &ypos);
		if (togglePicking) {
			pickingPos = { xpos, ypos };
			picking = true;
		}
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

int main() {

	#pragma region Setup

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
	const char *filename = "../models/scene.obj";
	Model model(filename);

	/////////////////////// Create cube and lines for cube
	float cubeSize = 2.f;
	Cube cube(model.center, cubeSize);
	GLuint cubeVAO = cube.vaoGeneration();
	std::cout << "Generated viewing cube." << std::endl;

	/////////////////////// Second cube for longer lines
	cubeSize = 20.f;
	Cube cube2(model.center, cubeSize);

	////////////////////// Create show bounding box
	GLuint bboxVAO = model.boundingCube.vaoGeneration();

	/////////////////////// Create show sphere
	Sphere sphere(model.center, model.radius);
	sphere.vaoGeneration(10, 20);
	Sphere sampleSphere(model.center, 2.f);
	sampleSphere.vaoGeneration(10, 20);

	///////////////////// Create sample locations
	int noSamples = 1000;
	int createRatio = 10;
	SphereSampler sampler(&sampleSphere, noSamples, createRatio);
	char dir = 'X';
	sampler.createSamplesFull();// createSamples(-1, dir); 
	sampler.vaoGeneration();
	glm::vec3 mainDir(1, 0, 0); /// check thisssss!!!!!!!!!!!
	std::cout << "Created " << sampler.samples.size() << " camera locations in " << dir << " direction" << std::endl;

	/////////////////// Create main camera
	PreviewCamera mainCamera;
	mainCamera.aspect = width / (float)height;

	////////////////// Create sampling canvas and sampling camera
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	GLuint rsttex, depthtex, framebuffer;
	const int w = 400, h = 400;
	glm::ivec2 resolution(w, h);
	setupTextureRender(rstProgram, w, h, model.vao, rsttex, depthtex, framebuffer);
	Orthocamera ocam(model.radius, model.radius, 0.1f, 30.f);
	Camera* cam = &ocam;

	#pragma endregion

	/////////////////// Make rayspacetree
	RaySpaceTree rst = RaySpaceTree();
	rst.depth = 8;
	rst.model = &model;
	int constructOption = 3;
	std::cout << "Constructing RST with depth = " << rst.depth << " and camera sample size = " << w << "x" << h << "..." << std::endl;
	//makeAdaptiveRST(&rst, sampler, rstProgram, resolution, cam, rsttex, constructOption);
	bool storeRays = true;
	makeRST(&rst, sampler, rstProgram, resolution, cam, rsttex, constructOption, storeRays);
	//fillMoreSamples(&rst, sampler, rsttex, rstProgram, resolution, cam, storeRays);

	////////////////// Statistics and checks
	getRstStatistics(&rst, model.primsize);
	//rst.checkLeaves();
	if (print) rst.printTree();

	////////////////// Check result by RayTracing
	if (trace) {
		glm::ivec2 traceRes = glm::ivec2(400, 400);
		cam->setPositionAndForward(glm::vec3(-2.f, 1.f, 0.f), model.center);
		checkTreeByTracing(cam, traceRes, &rst, rstProgram, rsttex);
	}

	#pragma region rendersetup

	////////////////// Cleanup
	glDeleteFramebuffers(1, &framebuffer);
	glDeleteTextures(1, &rsttex);
	glDeleteTextures(1, &depthtex);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	////////////////// Objects to intersect lines with
	GeoObject* geoObject;
	if (sphereOrCube) geoObject = &sphere;
	else geoObject = &cube;

	///////////////// Line Models & settings
	LineModel splitters = LineModel();
	splitters.updateVaoWithLines(rst.getAllSplittingLines(), geoObject, mainDir);
	LineModel samples = LineModel(true);
	LineModel extremalStabbing = LineModel();
	LineModel clickRay = LineModel();
	int leafnum = -1;
	Node* leaf;
	bool showing = false;
	bool viewline = false;
	float alphaLines = 1.f;

	////////////////// Colors
	glm::vec3 green = { 0, 1, 0 };
	glm::vec3 blue = { 0, 0, 1 };
	glm::vec3 red = { 1, 0, 0 };
	glm::vec3 yellow = { 1, 1, 0 };
	glm::vec3 white = { 1, 1, 1 };
	glm::vec3 black = { 0, 0, 0 };

	// Main loop
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	#pragma endregion

	while (!glfwWindowShouldClose(window)) {

		glfwPollEvents();

		if (nextleaf) {
			alphaLines = 1.f;
			int trinum = 0;
			while (trinum == 0) {
				leafnum++;
				leafnum = leafnum % rst.noLeaves;
				trinum = rst.getNumberOfTriInleaf(leafnum);
			}
			leaf = rst.getLeafFromNum(leafnum);
			nextleaf = false;
			update = true;
			viewline = false;
		}

		if (picking) {
			setupTextureRender(rstProgram, width, height, model.vao, rsttex, depthtex, framebuffer);

			GLfloat* pixels = cameraSamples(&mainCamera, rstProgram, rsttex, &model, glm::ivec2(width,height), 0);
			//int id = pixels
			glDeleteFramebuffers(1, &framebuffer);
			glDeleteTextures(1, &rsttex);
			glDeleteTextures(1, &depthtex);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			pickingPos = { 0, 0 };
			picking = false;
		}
		
		if (drawRay) {
			alphaLines = 0.5f;
			Ray r = mainCamera.pixRayDirection(drawRayPos);
			std::vector<Ray> rays = { r };
			leaf = rst.descend(r);
			//rst.descendWithLines(r, splitters, leafRays, rayColors, geoObject);
			glm::vec3 color(0, 1, 0);
			clickRay.updateVaoWithLines(rays, geoObject, mainDir);
			//rayVAO = Models::rayVao(r, color, geoObject, mainDir);
			drawRay = false;
			viewline = true;
			update = true;
		}

		if (update && leafnum >= 0) {
			if (sphereOrCube) geoObject = &sampleSphere;
			else geoObject = &cube;

			if (storeRays && toggleLines) 
				samples.updateVaoWithLines(rst.getViewingLinesInLeaf(leaf), geoObject, mainDir);
			if (toggle4lines)
				extremalStabbing.updateVaoWithLines(rst.getExtremalStabbingInLeaf(leaf), geoObject, mainDir);

			std::vector<Ray> split;
			rst.getSplittingLinesInLeaf(leaf, split);
			splitters.updateVaoWithLines(split, geoObject, mainDir);

			for (int i = 0; i < model.vertices.size(); i++) {
				model.vertices[i].selected = 0.f;
			}
			for (int i : leaf->primitiveSet) {
				model.vertices[3 * i].selected = 1.f;
				model.vertices[3 * i + 1].selected = 1.f;
				model.vertices[3 * i + 2].selected = 1.f;
			}
			model.changeSelected();
			update = false;
		}
		else update = false;

		updateCamera(mainCamera, width, height);
		glm::mat4 mvp = mainCamera.vpMatrix();

		// Bind the shader
		glUseProgram(mainProgram.index);
		glUniformMatrix4fv(glGetUniformLocation(mainProgram.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

		// Set viewport size
		glViewport(0, 0, width, height);

		// Clear the framebuffer to black and depth to maximum value
		glClearDepth(1.0f);  
		glClearColor(0.1f, 0.2f, 0.4f, 1.0f); 
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
		glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), 1.f);
		glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 0);

		if (toggleLines) {
			glBindVertexArray(splitters.vao);
			glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 1);
			glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(black));
			glDrawArrays(GL_LINES, 0, splitters.size);

			if (toggle4lines) {
				glBindVertexArray(extremalStabbing.vao);
				glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 1);
				glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(green));
				glDrawArrays(GL_LINES, 0, extremalStabbing.size);
				//glBindVertexArray(checklinesVao);
				//glm::vec3 blue = { 0, 1, 0 };
				//glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(blue));
				//glDrawArrays(GL_LINES, 0, through4lines.size());
				//glBindVertexArray(edgesVao);
				//glm::vec3 green = { 0, 1 , 0 };
				//glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(green));
				//glDrawArrays(GL_LINES, 0, 6);
			}
			if (toggleSampleLines) {
				glBindVertexArray(samples.vao);
				glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 0);
				glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), (GLfloat)alphaLines);
				glDrawArrays(GL_LINES, 0, samples.size);
			}
		}

		if (viewline) {
			glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 1);
			glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(red));
			glUniform1f(glGetUniformLocation(lineProgram.index, "alpha"), 1.f);
			glBindVertexArray(clickRay.vao);
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