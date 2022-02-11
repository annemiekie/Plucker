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
	bool alldir, trace, exact, sampling, cacheEE, cacheEEE, cacheCombi, rasterization, storeSamples, storeAllSamples;
	char dir;
	int sgn, depth, noSamples, createRatio, w, h, width, height;
	Options::constructOption construct;

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
		std::string constructStr = config.lookup("constructOption");
		if (Options::table.find(constructStr) != Options::table.end())
			construct = Options::table.find(constructStr)->second;
		else construct = Options::RANDOM_EDGE;
		std::string str = config.lookup("filename");
		filestr = str;
		width = config.lookup("width");
		height = config.lookup("height");
		trace = config.lookup("trace");
		sampling = config.lookup("sampling");
		exact = config.lookup("exact");
		cacheEE = config.lookup("cacheEE");
		cacheEEE = config.lookup("cacheEEE");
		cacheCombi = config.lookup("cacheCombi");
		rasterization = config.lookup("rasterization");
		storeSamples = config.lookup("storeSamples");
		storeAllSamples = config.lookup("storeAllSamples");

	}
	catch (libconfig::SettingNotFoundException& e) {
		std::cerr << "Incorrect setting(s) in configuration file." << std::endl;
	}

	Options::BuildOptions options = { construct, h, w, noSamples, rasterization, storeSamples, storeAllSamples, cacheCombi };

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
	if (sampling) rst = RSTBuilder<RSTBuilderSamples>::build(&model, depth, alldir, dir, sgn, options, visComp);
	model.enlargeModel();
	if (!sampling && exact)	rst = RSTBuilder<RSTBuilderExact>::build(&model, depth, alldir, dir, sgn, options, visComp);
	else if (exact) RSTBuilderExact::fill(&rst);

	rst.printLeafNodes();
	Visualizer::visualize(width, height, &model, &rst, visComp, window);


	return 0;
}