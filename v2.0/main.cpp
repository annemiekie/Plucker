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

#include <cstring>


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
#include "textureRenderer.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>

#include <CImg.h>
#include <libconfig.h++>

// Configuration
int width = 800;	// TESSST
int height = 800;	 // TESSST
bool update = false;
bool update4lines = false;
//bool trace = false;
bool togglePoints = false;
bool toggle4lines = false;
bool toggleSphere = false;
bool toggleBbox = false;
bool toggleLines = true;
bool toggleSampleLines = false;
bool nextleaf = false;
bool prevleaf = false;
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

void makeSample(glm::ivec2 res, Camera* cam, GLfloat* pixels, RaySpaceTree* rst, bool storeRays) {

	for (int y = 0; y < res.y; y++) {
		int yind = (res.y - 1 - y);
		for (int x = 0; x < res.x; x++) {
			const glm::vec2 pixpos { (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
			Ray ray = cam->pixRayDirection(pixpos);
			float tri = pixels[yind * res.x + x];
			if (tri > 0) {
				int tri_id = int(tri - 1);
				if (checkTriIntersect(tri_id * 3, ray, rst)) {
					rst->putPrimitive(ray, tri_id, storeRays);
				}
				//else if (tri_id2 >= 0) {
				else {
					//rst->putPrimitive(ray, tri_id2, storeRays);
					float t;
					if (rst->model->getIntersectionEmbree(ray, tri_id, t))
						rst->putPrimitive(ray, tri_id, storeRays);
				}

			}
		}
	}
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
			if (tri_id >= 0) {
				if (checkTriIntersect(tri_id * 3, ray, rst)) tris.insert(tri_id);
				else {
					float t;
					if (rst->model->getIntersectionEmbree(ray, tri_id, t))
						if (tri_id >= 0) tris.insert(tri_id);
				}
				samples.push_back({ count, tri_id});
			}
			else if (tri_id == -1) samples.push_back({ count, -1 });
			count++;
		}
	}
}

void fillMoreSamples(RaySpaceTree* rst, SphereSampler& sampler, Camera* cam, TextureRenderer& texrender) {

	std::cout << "Rasterizing camera samples to fill tree..." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < sampler.samples.size(); i++) {
		if (i % sampler.ratio == 0) continue;
		cam->setPositionAndForward(sampler.samples[i], rst->model->center);
		GLfloat* pixels = texrender.render(cam);
		makeSample(texrender.res, cam, pixels, rst, false);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	int noSamples = (sampler.samples.size() - (sampler.samples.size() / sampler.ratio + 1)) * texrender.height* texrender.width;
	std::cout << "Put " << noSamples << " samples in RST in " << diff << " ms." << std::endl;

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

void makeEmptyRST(RaySpaceTree* rst, int option) {
	std::vector<Ray> rays = constructRaysRandom(rst->model, rst->depth);
	rst->construct(option, rays);
}

void makeRST(RaySpaceTree* rst, SphereSampler& sampler, Camera* cam, TextureRenderer& texrender, int option, bool storeRays) {

	std::vector<Ray> rays = constructRaysRandom(rst->model, rst->depth);
	// Initialize tree
	rst->construct(option, rays);

	// Generate cameras
	std::vector<Ray> raytraceRays;

	std::cout << "Rasterizing camera samples to create tree..." << std::endl;

	for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
		cam->setPositionAndForward(sampler.samples[i], rst->model->center);
		GLfloat* pixels = texrender.render(cam);
		makeSample(texrender.res, cam, pixels, rst, storeRays);
	}
}

void makeRSTembree(RaySpaceTree* rst, SphereSampler& sampler, Camera* cam, int pixels_W, int pixels_H, int option, bool storeRays) {

	std::vector<Ray> rays = constructRaysRandom(rst->model, rst->depth);
	// Initialize tree
	rst->construct(option, rays);

	// Generate cameras
	std::vector<Ray> raytraceRays;

	std::cout << "Raytracing camera samples to create tree..." << std::endl;

	for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
		cam->setPositionAndForward(sampler.samples[i], rst->model->center);
		for (int y = 0; y < pixels_H; y++) {
			int yind = (pixels_H - 1 - y);
			for (int x = 0; x < pixels_W; x++) {
				const glm::vec2 pixpos{ (x + .5f) / pixels_W * 2.0f - 1.0f, 1.0f - (y + .5f) / pixels_H * 2.0f };
				Ray ray = cam->pixRayDirection(pixpos);
				int tri_id = -1;
				float t = 0.f;
				if (rst->model->getIntersectionEmbree(ray, tri_id, t))
					rst->putPrimitive(ray, tri_id, storeRays);
				//else rst->putPrimitive(ray, tri_id, storeRays, false);
			}
		}
	}
}

void makeAdaptiveRST(RaySpaceTree* rst, SphereSampler& sampler, Camera* cam, TextureRenderer& texrender, int option) {
	int imgsize = texrender.width * texrender.height;
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
		
		GLfloat* pixels = texrender.render(cam1);
		makeSample(texrender.res, cam1, pixels, rst, samples, tris, count);
		cams.push_back(cam);
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Made and stored " << samples.size() << " samples in " << diff << " ms." << std::endl;
	std::cout << "Removed " << texrender.width * texrender.height * (sampler.samples.size() / sampler.ratio + 1) - samples.size() << " samples for being too close to a primitive edge." << std::endl;

	std::cout << "Constructing the tree..." << std::endl;
	start_time = std::chrono::high_resolution_clock::now();
	if (option == 0) rst->constructAdaptive(texrender.res, cams, samples, tris, print);
	else if (option == 1) rst->constructSmartRandom(texrender.res, cams, samples, tris);
	end_time = std::chrono::high_resolution_clock::now();
	diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Constructed RST in " << diff << " ms." << std::endl;
}

void compareIntersectionMethods(Camera* cam, Model* model, TextureRenderer& texrender) {
	/////////////// compare with raytrace
	//GLfloat* pixels = texrender.render(cam);

	for (int y = 0; y < texrender.height; y++) {
		int yind = (texrender.height - 1 - y);
		for (int x = 0; x < texrender.width; x++) {
			const glm::vec2 pixpos{ (x + .5f) / texrender.width * 2.0f - 1.0f, 1.0f - (y + .5f) / texrender.width * 2.0f };
			Ray ray = cam->pixRayDirection(pixpos);
			//int triRast = int(pixels[yind * texrender.width + x])-1;
			float t1 = 0.f; float t2 = 0.f;
			int triTracePlucker = -1;
			model->getIntersectionNoAcceleration(ray, triTracePlucker, t1);
			int triTraceEmbree = -1;
			model->getIntersectionEmbree(ray, triTraceEmbree, t2);
			//if (triRast != triTracePlucker || triRast != triTraceEmbree || 
			if (triTraceEmbree != triTracePlucker) {
				//std::cout << "Rasterization: " << triRast << 
				std::cout << "Trace Plucker: " << triTracePlucker << " Trace Embree: " << triTraceEmbree << std::endl;
			}
			if (triTraceEmbree >= 0 && std::fabs(t1-t2) > 1E-8 ) {
				std::cout << "Trace Plucker depth: " << t1 << "Trace Embree depth: " << t2 << std::endl;
			}

		}
	}
}

void checkTreeByTracing(Camera* cam, RaySpaceTree* rst, TextureRenderer& texrender) {
	std::cout << "Raytracing to check tree..." << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();

	std::vector<int> raytraceTri = RayTracer::imageTracer(cam, rst, texrender.res);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Raytraced image in " << diff << " ms." << std::endl;

	/////////////// compare with raytrace
	GLfloat* pixels = texrender.render(cam);
	cimg_library::CImg<unsigned char> image(texrender.width, texrender.height, 1, 3, 0);
	const unsigned char white[] = { 255 , 255, 255 };
	const unsigned char black[] = { 0 , 0, 0 };

	for (int y = 0; y < texrender.height; y++) {
		int yind = (texrender.height - 1 - y);
		for (int x = 0; x < texrender.width; x++) {
			int triRast = int(pixels[yind * texrender.width + x]);
			int triTrace = raytraceTri[y * texrender.width + x];
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

////// MOVE TO RST
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
		case GLFW_KEY_0:
			prevleaf = true;
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
			//update4lines = true;
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
			drawRayPos = { (xpos + .5f) / float(width) * 2.0f - 1.0f, 1.0f - (ypos + .5f) / float(height) * 2.0f };
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
	/*create an instance of config*/
	libconfig::Config config;

	/*read a configuration file*/
	try {
		config.readFile("init.txt");
	}
	catch (libconfig::FileIOException& e) {
		/*inform user about IOException*/
		std::cerr << "FileIOException occurred. Could not read \"init.txt\"!!\n";
		/*terminate program*/
		exit(EXIT_FAILURE);
	}
	catch (libconfig::ParseException& e) {
		/*inform user about the parse exception*/
		std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine()
			<< " - " << e.getError() << std::endl;
		/*terminate program*/
		return(EXIT_FAILURE);
	}

	std::string filestr;// = "dragon.obj";
	bool alldir, trace, exact;
	char dir;
	int sgn, depth, noSamples, createRatio, w, h, constructOption;

	try {
		alldir = config.lookup("alldir");
		std::string dir1 = config.lookup("dir");
		dir = dir1.c_str()[0];
		sgn = config.lookup("sgn");
		depth = config.lookup("depth");
		noSamples = config.lookup("noSamples");
		createRatio = config.lookup("createRatio");
		w = config.lookup("w");
		h = config.lookup("h");
		constructOption = config.lookup("constructOption");
		std::string str = config.lookup("filename");
		filestr = str;
		width = config.lookup("width");
		height = config.lookup("height");
		trace = config.lookup("trace");
		exact = config.lookup("exact");
	}
	catch (libconfig::SettingNotFoundException& e) {
		std::cerr << "Incorrect setting(s) in configuration file." << std::endl;
	}


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

	GLFWwindow* window = glfwCreateWindow(width, height, "Ray Space Tree", nullptr, nullptr);
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
	const char* filename = filestr.c_str();
	Model model(filename, true);

	////////////////////// Create bounding box
	GLuint bboxVAO = model.boundingCube.vaoGeneration();

	/////////////////////// Create cube and lines for cube
	Cube cube(model.center, model.radius*1.5);
	GLuint cubeVAO = cube.vaoGeneration();
	std::cout << "Generated viewing cube." << std::endl;

	/////////////////////// Create show sphere
	Sphere sphere(model.center, model.radius);
	sphere.vaoGeneration(10, 20);
	Sphere sampleSphere(model.center, model.radius*1.5);
	sampleSphere.vaoGeneration(10, 20);

	///////////////////// Create sample locations
	SphereSampler sampler(&sampleSphere, noSamples, createRatio);
	if (alldir) dir = 'All';
	if (alldir) sampler.createSamplesFull();
	else sampler.createSamples(sgn, dir); 
	sampler.vaoGeneration();
	std::cout << "Created " << sampler.samples.size() << " camera locations in all directions" << std::endl;

	/////////////////// Create main camera
	PreviewCamera mainCamera(model.radius * 4, model.center, model.radius * 5 > 30.f ? model.radius * 5 : 30.f);
	mainCamera.aspect = width / (float)height;

	////////////////// Create sampling canvas and sampling camera
	TextureRenderer texrender = TextureRenderer(rstProgram, w, h, &model);
	Orthocamera ocam(model.radius, model.radius, model.radius * 2 > 30.f ? model.radius * 2 : 30.f);
	Camera* cam = &ocam;

	//cam->setPositionAndForward(glm::vec3(-2.f, 1.f, 0.f), model.center);
	//compareIntersectionMethods(cam, &model, texrender);

	#pragma endregion

	/////////////////// Make rayspacetree
	std::cout << "Constructing RST for direction " << dir << " with depth = " << depth << " and camera sample size = " << w << "x" << h << "..." << std::endl;
	RaySpaceTree rst = RaySpaceTree(&model, depth, alldir, sgn * glm::ivec3(dir == 'X', dir == 'Y', dir == 'Z'));
	bool storeRays = true;

	auto start_time = std::chrono::high_resolution_clock::now();
	
	if (!exact) makeRSTembree(&rst, sampler, cam, w, h, constructOption, storeRays);
	model.enlargeModel();
	if (exact) {
		makeEmptyRST(&rst, constructOption);
		rst.fillExact();
	}

	for (Node* n : rst.nodes) {
		if (n->leaf) {
			std::cout << n->index-rst.noLeaves+1 << ": ";
			for (int t : n->primitiveSet)
				std::cout << t << " ";
			std::cout << std::endl;
		}
	}
	
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Completed RST1 in " << diff << " ms" << std::endl;

	////////////////// Statistics and checks
	getRstStatistics(&rst, model.primsize);

	if (print) rst.printTree();

	////////////////// Check result by RayTracing
	if (trace) {
		glm::ivec2 traceRes = glm::ivec2(400, 400);
		cam->setPositionAndForward(glm::vec3(-2.f, 1.f, 0.f), model.center);
		checkTreeByTracing(cam, &rst, texrender);
	}

	#pragma region rendersetup

	////////////////// Cleanup
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	////////////////// Objects to intersect lines with
	GeoObject* geoObject;
	if (sphereOrCube) geoObject = &sphere;
	else geoObject = &cube;

	///////////////// Line Models & settings
	LineModel splitters = LineModel();
	splitters.updateVaoWithLines(rst.getAllSplittingLines(), geoObject, rst.maindir);
	LineModel samples = LineModel(true);
	LineModel extremalStabbing = LineModel();
	LineModel clickRay = LineModel();
	LineModel wrongRays = LineModel();
	int leafnum = -1;
	Node* leaf;
	bool showing = false;
	bool viewline = false;
	float alphaLines = 1.f;

	GLuint sideQuadVao = 0;
	if (!alldir) sideQuadVao = model.boundingCube.vaoSideQuad(rst.maindir);

	////////////////// Colors
	glm::vec3 green = { 0, 1, 0 };
	glm::vec3 blue = { 0, 0, 1 };
	glm::vec3 red = { 1, 0, 0 };
	glm::vec3 yellow = { 1, 1, 0 };
	glm::vec3 white = { 1, 1, 1 };
	glm::vec3 black = { 0, 0, 0 };
	glm::vec3 hotpink = { 1, 0.4, 0.7 };

	// Main loop
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Lighting
	float lightdis = 3 * model.radius;
	glm::vec3 lightPos = glm::vec3(-lightdis, 2*lightdis, -lightdis);
	glm::vec3 lightPos2 = glm::vec3(-lightdis, 2*lightdis, lightdis);
	glm::vec3 lightPos3 = glm::vec3(lightdis, -2*lightdis, 0);
	glUseProgram(mainProgram.index);
	glUniform3fv(glGetUniformLocation(mainProgram.index, "lightPos"), 1, glm::value_ptr(lightPos));
	glUniform3fv(glGetUniformLocation(mainProgram.index, "lightPos2"), 1, glm::value_ptr(lightPos2));
	glUniform3fv(glGetUniformLocation(mainProgram.index, "lightPos3"), 1, glm::value_ptr(lightPos3));
	#pragma endregion


	//rst.checkLeaves();


	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		glUseProgram(0);
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

		if (prevleaf) {
			alphaLines = 1.f;
			int trinum = 0;
			while (trinum == 0) {
				leafnum--;
				leafnum = leafnum % rst.noLeaves;
				trinum = rst.getNumberOfTriInleaf(leafnum);
			}
			leaf = rst.getLeafFromNum(leafnum);
			prevleaf = false;
			update = true;
			viewline = false;
		}

		if (picking) {
			Ray r = mainCamera.pixRayDirection(drawRayPos);
			int primindex = -1;
			float t = 0;
			//model.getIntersectionNoAcceleration(r, primindex, t);
			model.getIntersectionEmbree(r, primindex, t);
			if (primindex >= 0) {
				std::cout << primindex << std::endl;
				Ray extremalLine;
				if (rst.check1Prim(primindex, extremalLine, leaf, true, 0))
					extremalStabbing.updateVaoWithLines(std::vector<Ray>{extremalLine}, geoObject, rst.maindir);
			}
			//model.vertices[3 * primindex].selected = 1.f;
			//model.vertices[3 * primindex + 1].selected = 1.f;
			//model.vertices[3 * primindex + 2].selected = 1.f;
			picking = false;
		}
		
		if (drawRay) {
			alphaLines = 0.5f;
			Ray r = mainCamera.pixRayDirection(drawRayPos);
			std::vector<Ray> rays = { r };
			leaf = rst.descend(r);
			glm::vec3 color(0, 1, 0);
			clickRay.updateVaoWithLines(rays, geoObject, rst.maindir);
			drawRay = false;
			viewline = true;
			update = true;
		}

		if (update && leafnum >= 0) {
			std::cout << leafnum << std::endl;

			if (sphereOrCube) geoObject = &sampleSphere;
			else geoObject = &cube;

			std::vector<Ray> split;
			rst.getSplittingLinesInLeaf(leaf, split);
			splitters.updateVaoWithLines(split, geoObject, rst.maindir);

			for (int i = 0; i < model.vertices.size(); i++) {
				model.vertices[i].selected = 0.f;
			}
			for (int i : leaf->primitiveSet) {
				model.vertices[3 * i].selected = 1.f;
				model.vertices[3 * i + 1].selected = 1.f;
				model.vertices[3 * i + 2].selected = 1.f;
			}
			model.changeSelected();

			if (storeRays && toggleLines)
				samples.updateVaoWithLines(rst.getViewingLinesInLeaf(leaf), geoObject, rst.maindir);
			if (toggle4lines) {// && update4lines) {
				extremalStabbing.updateVaoWithLines(rst.getExtremalStabbingInLeaf(leaf, true), geoObject, rst.maindir);
				//update4lines = false;
				//wrongRays.updateVaoWithLines(rst.wronglines, geoObject, rst.maindir);
			}


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
			if (!alldir) {
				glBindVertexArray(sideQuadVao);
				glDrawArrays(GL_LINES, 0, 8);
			}
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
			glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(yellow));
			glDrawArrays(GL_LINES, 0, splitters.size);

			if (toggle4lines) {
				//glBindVertexArray(extremalStabbing.vao);
				//glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 1);
				//glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(green));
				//glDrawArrays(GL_LINES, 0, extremalStabbing.size);

				glBindVertexArray(extremalStabbing.vao);
				glUniform1i(glGetUniformLocation(lineProgram.index, "setcol"), 1);
				glUniform3fv(glGetUniformLocation(lineProgram.index, "setcolor"), 1, glm::value_ptr(green));
				glDrawArrays(GL_LINES, 0, extremalStabbing.size);

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