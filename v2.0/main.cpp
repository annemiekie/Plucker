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

#include "rst.h"
#include "model.h"
#include "rstBuilder.h"
#include "rstBuilderSamples.h"
#include "rstBuilderExact.h"
#include "buildOptions.h"
#include "visComponents.h"
#include "visualizer.h"
#include "axisAlignedPolygon.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>

#include <CImg.h>
#include <libconfig.h++>

// OpenGL debug callback
void APIENTRY debugCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {
	if (severity != GL_DEBUG_SEVERITY_NOTIFICATION) {
		std::cerr << "OpenGL: " << message << std::endl;
	}
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

	std::string filestr, constructStr;
	char dir;
	int sgn, width, height;
	Options::BuildOptions options;

	try {
		options.alldir = config.lookup("alldir");
		std::string dir1 = config.lookup("dir");
		dir = dir1.c_str()[0];
		sgn = config.lookup("sgn");
		options.depth = config.lookup("depth");

		std::string str = config.lookup("filename");
		filestr = str;
		width = config.lookup("width");
		height = config.lookup("height");

		std::string constructStr = config.lookup("constructOption");
		if (Options::table.find(constructStr) != Options::table.end())
			options.construct = Options::table.find(constructStr)->second;
		else options.construct = Options::RANDOM_EDGE;

		std::string samplingStr = config.lookup("samplingType");
		if (Options::table2.find(samplingStr) != Options::table2.end())
			options.samplingtype = Options::table2.find(samplingStr)->second;
		else options.samplingtype = Options::UNIFORM_SPHERE;

		options.sampling = config.lookup("sampling");
		options.noSamples = config.lookup("noSamples");
		options.sampleStoreFillRatio = config.lookup("sampleStoreFillRatio");
		options.s_w = config.lookup("w");
		options.s_h = config.lookup("h");
		options.rasterizationSampling = config.lookup("rasterization");
		options.storeSamples = config.lookup("storeSamples");
		options.storeAllSamples = config.lookup("storeAllSamples");
		options.fillMoreSamples = config.lookup("fillMoreSamples");

		options.exact = config.lookup("exact");
		options.exactStartLevel = config.lookup("exactStartLevel");
		options.exactStartLevel = std::min(options.depth, options.exactStartLevel);

		options.cacheEE = config.lookup("cacheEE");
		options.cacheEEE = config.lookup("cacheEEE");
		options.cacheCombi = config.lookup("cacheCombi");


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

	// Activate the OpenGL context
	glfwMakeContextCurrent(window);

	// Initialize GLEW extension loader
	glewExperimental = GL_TRUE;
	glewInit();

	// Set up OpenGL debug callback
	glDebugMessageCallback(debugCallback, nullptr);

	////////////////////////// Load vertices of model
	const char* filename = filestr.c_str();
	Model model(filename);

	VisComponents visComp;
	RaySpaceTree rst;
	glm::vec3 maindir = glm::vec3(0);
	if (!options.alldir) maindir = sgn * glm::ivec3(dir == 'X', dir == 'Y', dir == 'Z');
	AxisAlignedPolygon polyWindow = model.boundingCube.getCubeSideSquare(maindir);
	//polyWindow.makeSmaller(0.2);

	if (options.sampling) rst = RSTBuilder<RSTBuilderSamples>::build(&model, &polyWindow, dir, sgn, options, visComp, false);
	else if (options.exact)	rst = RSTBuilder<RSTBuilderExact>::build(&model, &polyWindow, dir, sgn, options, visComp, true);
	if (options.sampling && options.exact) RSTBuilderExact::fill(options, &rst, true);
	//rst.printLeafNodes();
	ESLFindAll eslfinder(&rst);
	eslfinder.find();

	Visualizer::visualize(width, height, &model, &rst, visComp, window, eslfinder.raysPerPrimitive);


	return 0;
}