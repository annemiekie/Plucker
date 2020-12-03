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

// Library for loading .OBJ model
#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>

// Library for loading an image
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>


// Header for camera structure/functions
#include "camera.h"
#include "ray.h"
#include "vertex.h"
#include "rst.h"
#include "cube.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <CImg.h>


// Configuration
const int width = 800;
const int height = 800;
int leafNumber = -1;

// Helper function to read a file like a shader
std::string readFile(const std::string& path) {
	std::ifstream file(path, std::ios::binary);
	
	std::stringstream buffer;
	buffer << file.rdbuf();

	return buffer.str();
}

bool checkShaderErrors(GLuint shader) {
	// Check if the shader compiled successfully
	GLint compileSuccessful;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compileSuccessful);

	// If it didn't, then read and print the compile log
	if (!compileSuccessful) {
		GLint logLength;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);

		std::vector<GLchar> logBuffer(logLength);
		glGetShaderInfoLog(shader, logLength, nullptr, logBuffer.data());

		std::cerr << logBuffer.data() << std::endl;
		
		return false;
	} else {
		return true;
	}
}

bool checkProgramErrors(GLuint program) {
	// Check if the program linked successfully
	GLint linkSuccessful;
	glGetProgramiv(program, GL_LINK_STATUS, &linkSuccessful);

	// If it didn't, then read and print the link log
	if (!linkSuccessful) {
		GLint logLength;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);

		std::vector<GLchar> logBuffer(logLength);
		glGetProgramInfoLog(program, logLength, nullptr, logBuffer.data());

		std::cerr << logBuffer.data() << std::endl;
		
		return false;
	} else {
		return true;
	}
}

// OpenGL debug callback
void APIENTRY debugCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {
	if (severity != GL_DEBUG_SEVERITY_NOTIFICATION) {
		std::cerr << "OpenGL: " << message << std::endl;
	}
}

GLuint generateShader(std::string vert, std::string frag) {
	////////////////// Load and compile rst shader program
	GLuint shaderProgram = glCreateProgram();

	std::string vertexShaderCode = readFile(vert);
	const char* vertexShaderCodePtr = vertexShaderCode.data();

	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexShaderCodePtr, nullptr);
	glCompileShader(vertexShader);

	std::string fragmentShaderCode = readFile(frag);
	const char* fragmentShaderCodePtr = fragmentShaderCode.data();

	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentShaderCodePtr, nullptr);
	glCompileShader(fragmentShader);

	if (!checkShaderErrors(vertexShader) || !checkShaderErrors(fragmentShader)) {
		std::cerr << "Shader(s) failed to compile!" << std::endl;
		return EXIT_FAILURE;
	}

	// Combine vertex and fragment shaders into single shader program
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);

	if (!checkProgramErrors(shaderProgram)) {
		std::cerr << "Shadow program failed to link!" << std::endl;
		return EXIT_FAILURE;
	}
	return shaderProgram;
}

// Key handle function
void keyboardHandler(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	cameraKeyboardHandler(key, action);
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


GLuint loadModel(std::vector<Vertex>& vertices, int& size, const char* filename) {

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;

	if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename)) {
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	size = shapes[1].mesh.indices.size();
	int ind = 0;
	int count = 0;
	// Read triangle vertices from OBJ file
	//for (const auto& shape : shapes) {
	for (const auto& index : shapes[1].mesh.indices) {
		Vertex vertex = {};

		// Retrieve coordinates for vertex by index
		vertex.pos = {
			attrib.vertices[3 * index.vertex_index + 0],
			attrib.vertices[3 * index.vertex_index + 1],
			attrib.vertices[3 * index.vertex_index + 2]
		};

		// Retrieve components of normal by index
		vertex.normal = {
			attrib.normals[3 * index.normal_index + 0],
			attrib.normals[3 * index.normal_index + 1],
			attrib.normals[3 * index.normal_index + 2]
		};

		if (count % 3 == 0) ind++;
		vertex.id = (1.f*ind) / size;
		count++;

		vertices.push_back(vertex);
	}
	//}

	//////////////////// Create Vertex Buffer Object
	GLuint vbo;
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

	// Bind vertex data to shader inputs using their index (location)
	// These bindings are stored in the Vertex Array Object
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
	// Stride is the distance in bytes between vertices
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, pos)));
	glEnableVertexAttribArray(0);

	// The normals should be retrieved from the same Vertex Buffer Object (glBindBuffer is optional)
	// The offset is different and the data should go to input 1 instead of 0
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, normal)));
	glEnableVertexAttribArray(1);

	// The normals should be retrieved from the same Vertex Buffer Object (glBindBuffer is optional)
	// The offset is different and the data should go to input 1 instead of 0
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, id)));
	glEnableVertexAttribArray(2);

	return vao;
}

GLuint generateDepthBuffer(const int rstwidth, const int rstheight) {
	GLuint depthrenderbuffer;
	glGenRenderbuffers(1, &depthrenderbuffer);
	glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, rstwidth, rstheight);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);
	return depthrenderbuffer;
}

GLuint generateFrameBuffer(GLuint texRST) {
	GLuint framebuffer;
	glGenFramebuffers(1, &framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	// Set "renderedTexture" as our colour attachement #0
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texRST, 0);
	return framebuffer;
}

GLuint generateTexture(const int rstwidth, const int  rstheight) {
	GLuint texRST;

	glGenTextures(1, &texRST);
	glBindTexture(GL_TEXTURE_2D, texRST);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, rstwidth, rstheight, 0, GL_RGB, GL_FLOAT, nullptr);

	//// Set interpolation for texture sampling (GL_NEAREST for no interpolation)
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);
	return texRST;
}

GLuint quadGeneration() {
	float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
	// positions   // texCoords
		-1.0f,  1.0f,  0.0f, 1.0f,
		-1.0f, -1.0f,  0.0f, 0.0f,
		1.0f, -1.0f,  1.0f, 0.0f,

		-1.0f,  1.0f,  0.0f, 1.0f,
		1.0f, -1.0f,  1.0f, 0.0f,
		1.0f,  1.0f,  1.0f, 1.0f
	};

	// quaaad
	unsigned int quadVAO, quadVBO;
	glGenVertexArrays(1, &quadVAO);
	glGenBuffers(1, &quadVBO);
	glBindVertexArray(quadVAO);
	glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));

	return quadVAO;
}

GLuint lineGeneration(glm::vec3& min, glm::vec3& max) {
	std::vector<GLfloat> lineSegments = {
		min.x, min.y, min.z,
		max.x, min.y, min.z,
		min.x, min.y, min.z,
		min.x, min.y, max.z,
		max.x, min.y, min.z,
		max.x, min.y, max.z,
		min.x, min.y, max.z,
		max.x, min.y, max.z,

		min.x, max.y, min.z,
		max.x, max.y, min.z,
		min.x, max.y, min.z,
		min.x, max.y, max.z,
		max.x, max.y, min.z,
		max.x, max.y, max.z,
		min.x, max.y, max.z,
		max.x, max.y, max.z,

		min.x, min.y, min.z,
		min.x, max.y, min.z,
		min.x, min.y, max.z,
		min.x, max.y, max.z,
		max.x, min.y, max.z,
		max.x, max.y, max.z,
		max.x, min.y, min.z,
		max.x, max.y, min.z,
	};

	//GLfloat lineSeg[] =
	//{
	//	-1.f, 0.f,-1.f,
	//	1.f, 0.f,-1.f,
	//	-1.f, 0.f, -1.f,
	//	-1.f, 0.f, 1.f,
	//	1.f, 0.f,1.f,
	//	1.f, 0.f,-1.f,
	//	1.f, 0.f,1.f,
	//	-1.f, 0.f,1.f,
//
	//	-1.f, 2.f,-1.f,
	//	1.f, 2.f,-1.f,
	//	-1.f, 2.f, -1.f,
	//	-1.f, 2.f, 1.f,
	//	1.f, 2.f,1.f,
	//	1.f, 2.f,-1.f,
	//	1.f, 2.f,1.f,
	//	-1.f, 2.f,1.f,
//
	//	-1.f, 0.f,1.f,
	//	-1.f, 2.f,1.f,
	//	1.f, 0.f,-1.f,
	//	1.f, 2.f,-1.f,
	//	-1.f, 0.f,-1.f,
	//	-1.f, 2.f,-1.f,
	//	1.f, 0.f,1.f,
	//	1.f, 2.f,1.f,
	//};
	//lines = 24;

	//for (int i = 0; i < lines*3; i++) {
	//	lineSeg[i] /= 2.f;
	//}

	GLuint lineVAO, lineVBO;
	glGenVertexArrays(1, &lineVAO);
	glGenBuffers(1, &lineVBO);
	glBindVertexArray(lineVAO);
	glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * lineSegments.size(), &lineSegments[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
	return lineVAO;
}

GLuint vaoLineGeneration(std::vector<glm::vec3> &lines) {

	GLuint lineVBO;
	glGenBuffers(1, &lineVBO);
	glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * lines.size(), &lines[0].x, GL_STATIC_DRAW);

	GLuint lineVAO;
	glGenVertexArrays(1, &lineVAO);
	glBindVertexArray(lineVAO);

	glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
	glEnableVertexAttribArray(0);

	return lineVAO;
}

GLuint vaoLineGenerationWithColor(std::vector<glm::vec3>& lines, std::vector<glm::vec3>& colorpos, Cube *cube) {

	std::vector<glm::vec3> colors;
	glm::vec3 color;
	glm::vec3 min = cube->bounds[0] - glm::vec3(cube->size);
	float size = cube->size * 3.f;

	for (int i = 0; i < lines.size() / 2; i++) {
		glm::vec3 pos = colorpos[i];
		color = (pos - min) / size;
		color.x = 0;
		
		colors.push_back(color);
		colors.push_back(color);
	}

	GLuint lineVBO[2];
	glGenBuffers(2, lineVBO);
	glBindBuffer(GL_ARRAY_BUFFER, lineVBO[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * lines.size(), &lines[0].x, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, lineVBO[1]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * colors.size(), &colors[0].x, GL_STATIC_DRAW);

	GLuint lineVAO;
	glGenVertexArrays(1, &lineVAO);
	glBindVertexArray(lineVAO);

	glBindBuffer(GL_ARRAY_BUFFER, lineVBO[0]);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, lineVBO[1]);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
	glEnableVertexAttribArray(1);

	return lineVAO;
}

// do this in the shader as well? Or cuda?
Ray raydirection(glm::ivec2 resolution, Camera& cam, glm::ivec2 pixel, glm::mat4& camToWorld, float halfScreenPlace) {
	// NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
	const glm::vec2 normalizedPixelPos{
		(pixel.x + .5f) / resolution.x * 2.0f - 1.0f,
		1.0f - (pixel.y + .5f) / resolution.y * 2.0f
	};
	glm::vec3 camSpaceRayDirection{
		normalizedPixelPos.x * halfScreenPlace,
		normalizedPixelPos.y * halfScreenPlace * cam.aspect,
		-1
	};
	camSpaceRayDirection = glm::normalize(camSpaceRayDirection);
	glm::vec3 worldSpaceRayDirection = glm::normalize(camToWorld * glm::vec4(camSpaceRayDirection, 0.0f));
	return Ray(cam.position + worldSpaceRayDirection, cam.position);
}

void makeSample(glm::ivec2 resolution, Camera& cam, std::vector<float>& pixels, 
				int triNum, RaySpaceTree* rst) {
	const float halfScreenPlace = std::tan(cam.fov / 2.0f);
	//should it also be projection matrix?
	glm::mat4 camToWorld = glm::inverse(cam.vMatrix());

	cimg_library::CImg<unsigned char> image(resolution.x, resolution.y, 1, 3, 0);
	for (int y = 0; y < resolution.y; y++) {
		for (int x = 0; x < resolution.x; x++) {
			Ray ray = raydirection(resolution, cam, glm::ivec2(x, y), camToWorld, halfScreenPlace);
			float tri = pixels[(resolution.y - 1 - y) * resolution.x + x];
			int tri_id = 0;
			if (tri > 0) {
				// always add ray to tree to calculate size of a node
				tri_id = int(triNum * tri) - 1;
				rst->putPrimitive(ray, tri_id);
				const unsigned char color[] = { 255 , 0, 0 };
				image.draw_point(x, y, color);
			}

		}
	}
	image.save("file.bmp");
}

void cameraSamples(Camera& cam, GLuint rstshader, GLuint rsttex, int no_triangles, 
					glm::ivec2 resolution, std::vector<Vertex>& vertices, RaySpaceTree* rst) {
	glm::mat4 mvp = cam.vpMatrix();
	glUniformMatrix4fv(glGetUniformLocation(rstshader, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
	glUniform3fv(glGetUniformLocation(rstshader, "viewPos"), 1, glm::value_ptr(cam.position));

	glClearDepth(1.0f);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	// Execute draw command
	glDrawArrays(GL_TRIANGLES, 0, vertices.size());

	int size = resolution.x * resolution.y;
	glBindTexture(GL_TEXTURE_2D, rsttex);
	GLfloat* pixels = new GLfloat[size * 3];
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pixels);
	std::vector<float> pix(size);
	for (int i = 0; i < size; i++)
		pix[i] = (float)pixels[i * 3];
	makeSample(resolution, cam, pix, no_triangles, rst);

}

void makeRST(RaySpaceTree* rst, std::vector<Vertex> &vertices, GLuint vao, int no_triangles, GLFWwindow* window, 
			GLuint shader, GLuint rstshader, int rstwidth, int rstheight, int option) {

	GLuint rsttex = generateTexture(rstwidth, rstheight);
	GLuint framebuffer = generateFrameBuffer(rsttex);
	GLuint depthbuffer = generateDepthBuffer(rstwidth, rstheight);
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) return;

	// Clear the map and set needed options
	glClearDepth(1.0f);
	glClear(GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	// Bind the shader
	glUseProgram(rstshader);
	glUniform1f(glGetUniformLocation(rstshader, "noPrim"), static_cast<float>(no_triangles));

	// Bind vertex data & buffer
	glBindVertexArray(vao);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

	// Set viewport size
	glViewport(0, 0, rstwidth, rstheight);

	// Initialize tree
	rst->depth = 8;
	rst->mainDir = glm::vec3(1, 0, 0);
	rst->construct(0, vertices, option);

	// Generate cameras
	Camera cam;
	cam.aspect = rstwidth/rstheight;

	int noSample = 10;
	float sampleWidth = 2.f / (float)noSample;
	float sampleStart[2] = { 0.f, -1.f };
	//need cube to check intersection - if no intersection, not in rst!!
	// Sample per cameras
	for (int z = 0; z < noSample; z++) {
		for (int y = 0; y < noSample; y++) {
			cam.position = glm::vec3(-2.f, sampleStart[0] + sampleWidth*y, sampleStart[0] + sampleWidth * z);
			cam.forward = -cam.position;
			cameraSamples(cam, rstshader, rsttex, no_triangles, glm::ivec2(rstwidth, rstheight), vertices, rst);
		}
	}

	//GLuint quadvao = quadGeneration();
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
		std::cout <<  "Press enter to close."; getchar();
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
	GLuint mainProgram = generateShader("standard.vert", "standard.frag");
	GLuint rstProgram = generateShader("rst.vert", "rst.frag");
	GLuint texprojProgram = generateShader("texproj.vert", "texproj.frag");
	GLuint lineProgram = generateShader("line.vert", "line.frag");

	////////////////////////// Load vertices of model
	std::vector<Vertex> vertices;
	int noPrimitives = 0;
	const char *filename = "scene.obj";
	glm::ivec4 includes( 12299, 35876, 6149, 27675 );
	GLuint vao = loadModel(vertices, noPrimitives, filename);

	/////////////////////// Create cube and lines for cube
	float cubeSize = 2.f;
	glm::vec3 cubeMin = glm::vec3(-1.f, 0.f, -1.f);
	glm::vec3 cubeMax = cubeMin + cubeSize;
	Cube cube(cubeMin, cubeSize);
	GLuint lineVAO = lineGeneration(cubeMin, cubeMax);

	/////////////////// Create main camera
	Camera mainCamera;
	mainCamera.aspect = width / (float)height;
	mainCamera.position = glm::vec3(1.1f, 0.9f, 0.9f);
	mainCamera.forward = glm::vec3(-1.1f, -0.5f, -0.9f);

	////////////////// Unbind the off-screen framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	/////////////////// Make rayspacetree
	std::vector<glm::vec3> splitters;
	std::vector<glm::vec3> leafRays;
	std::vector<glm::vec3> rayColors;

	const int w = 400, h = 400;
	RaySpaceTree rst(&cube);
	makeRST(&rst, vertices, vao, noPrimitives, window, texprojProgram,
								rstProgram, w, h, 0);
	// Use rst info to visualize
	splitters = rst.getSplittingLinesInCube();

	GLuint splitterVao = vaoLineGeneration(splitters);
	GLuint leafRayVao;
	int noLeaves = std::pow(2, rst.depth);
	int leafnum = leafNumber;
	bool showing = false;

	//getRstStatistics(&rst, noPrimitives);

	// Main loop
	while (!glfwWindowShouldClose(window)) {

		glfwPollEvents();

		if (leafnum != leafNumber) {
			leafRays = std::vector<glm::vec3>();
			splitters = std::vector<glm::vec3>();
			rayColors = std::vector<glm::vec3>();

			while (leafRays.size() == 0) {
				leafNumber = leafNumber % noLeaves;
				rst.getViewingLinesInLeaf(leafNumber, leafRays, rayColors, splitters);
				splitterVao = vaoLineGeneration(splitters);
				leafNumber++;
			}

			//if (leafRays.size() > 0) {
				leafRayVao = vaoLineGenerationWithColor(leafRays, rayColors, &cube);
				showing = true;
			//}
			leafnum = leafNumber;
		}

		updateCamera(mainCamera, width, height);
		glm::mat4 mvp = mainCamera.vpmMatrix();

		// Bind the shader
		glUseProgram(mainProgram);
		glUniformMatrix4fv(glGetUniformLocation(mainProgram, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
		glUniform3fv(glGetUniformLocation(mainProgram, "viewPos"), 1, glm::value_ptr(mainCamera.position));
		glUniform1i(glGetUniformLocation(mainProgram, "noPrim"), noPrimitives);
		glUniform4iv(glGetUniformLocation(mainProgram, "includes"), 1, glm::value_ptr(includes));

		// Bind vertex data
		glBindVertexArray(vao);
	
		// Set viewport size
		glViewport(0, 0, width, height);

		// Clear the framebuffer to black and depth to maximum value
		glClearDepth(1.0f);  
		glClearColor(0.1f, 0.2f, 0.3f, 1.0f); 
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDisable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);

		// Execute draw command
		glDrawArrays(GL_TRIANGLES, 0, vertices.size());

		glBindVertexArray(lineVAO);
		glDrawArrays(GL_LINES, 0, 24);

		glBindVertexArray(splitterVao);
		glDrawArrays(GL_LINES, 0, splitters.size());

		if (showing) {
			glUseProgram(lineProgram);
			glUniformMatrix4fv(glGetUniformLocation(lineProgram, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			glBindVertexArray(leafRayVao);
			glDrawArrays(GL_LINES, 0, leafRays.size());
		}

		// Present result to the screen
		glfwSwapBuffers(window);
	}
	
	//glDeleteFramebuffers(1, &framebuffer);

	//glDeleteTextures(1, &texRST);

	glfwDestroyWindow(window);
	
	glfwTerminate();

    return 0;
}