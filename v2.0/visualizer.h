#pragma once
#define GLEW_STATIC
#include <GL/glew.h>

// Library for window creation and event handling
#include <GLFW/glfw3.h>

// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <vector>

#include "previewCamera.h"
#include "shader.h"
#include "model.h"
#include "linemodel.h"
#include "geoObject.h"
#include "rst.h"
#include "colors.h"
#include "lights.h"
#include "visComponents.h"

namespace Visualizer {
	int width, height;
	bool togglePoints, toggle4lines, toggleSphere, toggleBbox, toggleLines, toggleSampleLines, togglePicking, toggleGeo;
	GLFWwindow* window;
	Shader mainShader, lineShader;
	Model* model;
	PreviewCamera cam;
	LineModel splitters, samples, clickRay, edgeRays, eslSilhEdges, eslLines;
	GeoObject* geoObject;
	RaySpaceTree* rst;
	Node* leaf;
	int leafnum;
	bool showing, changeColoring, viewline, edgelines;
	float alphaLines;
	int selectedPrim;
	std::vector<int> notfoundprim;
	Lights lights;
	VisComponents visComp;
	//static Visualizer vis;


	void drawLines(LineModel& lines, int setcol, float alpha, glm::vec3 color = glm::vec3()) {
		glUniform1i(glGetUniformLocation(lineShader.index, "setcol"), setcol);
		if (setcol) glUniform3fv(glGetUniformLocation(lineShader.index, "setcolor"), 1, glm::value_ptr(color));
		glUniform1f(glGetUniformLocation(lineShader.index, "alpha"), alpha);
		lines.draw();
	}

	void drawElements() {

		glViewport(0, 0, width, height);
		glClearDepth(1.0f);
		glClearColor(0.5f, 0.6f, 0.65f, 1);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		glm::mat4 mvp = cam.vpMatrix();
		glUseProgram(mainShader.index);
		glUniformMatrix4fv(glGetUniformLocation(mainShader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
		model->draw();
		model->boundingCube.draw();

		if (togglePoints) {
			glBindVertexArray(visComp.sampling_vao);
			glDrawArrays(GL_POINTS, 0, visComp.sampling_size);
		}
		if (toggleBbox) model->boundingBox.draw();


		glUseProgram(lineShader.index);
		glUniformMatrix4fv(glGetUniformLocation(lineShader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
		if (toggleLines) {
			if (rst->splitters.size() > 0) drawLines(splitters, 1, 1.f, Colors::yellow);
			if (toggle4lines) drawLines(eslLines, 1, 1.f, Colors::hotpink);
			if (toggleSampleLines) drawLines(samples, 0, (GLfloat)alphaLines);
			if (viewline) drawLines(clickRay, 1, 1.f, Colors::red);
			if (edgelines) {
				drawLines(edgeRays, 1, 1.f, Colors::green);
				drawLines(eslSilhEdges, 1, 1.f, Colors::orange);
			}
		}

		if (toggleSphere) {
			mvp = cam.vpmMatrix(model->boundingSphere.center);
			glUseProgram(mainShader.index);
			glUniformMatrix4fv(glGetUniformLocation(mainShader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			model->boundingSphere.draw();
		}

	}

	void checkColorChange() {
		if (changeColoring) {
			model->clearColors(Colors::white);
			model->changePrimColors<std::set<int>>(leaf->primitiveSet, Colors::bluegreen);
			model->changePrimColors<std::vector<int>>(notfoundprim, Colors::red);
			if (selectedPrim >= 0) model->changePrimColor(selectedPrim, Colors::green);
			model->changeSelected();
			changeColoring = false;
		}
	}


	glm::vec2 getDrawRayPos(int xpos, int ypos) {
		return { (xpos + .5f) / float(width) * 2.0f - 1.0f, 1.0f - (ypos + .5f) / float(height) * 2.0f };
	}


	void updateLines() {
		if (leafnum < 0) return;

		std::vector<Ray> split;
		rst->getSplittingLinesInLeaf(leaf, split);
		splitters.updateVaoWithLines(split, geoObject);

		changeColoring = true;
		notfoundprim = std::vector<int>();
		if (toggleSampleLines)
			samples.updateVaoWithLines(rst->getViewingLinesInLeaf(leaf), geoObject, rst->maindir);
		if (toggle4lines)
			eslLines.updateVaoWithLines(rst->getExtremalStabbingInLeaf(leaf, notfoundprim), geoObject);
	}

	void leafChange(int num) {
		selectedPrim = -1;
		alphaLines = 1.f;
		int trinum = 0;
		while (trinum == 0) {
			leafnum += num;
			leafnum = leafnum % rst->noLeaves;
			trinum = rst->getNumberOfTriInleaf(leafnum);
		}
		leaf = rst->getLeafFromNum(leafnum);
		updateLines();
		viewline = false;
		edgelines = false;
	}

	void changeGeoObject(bool geosphere) {
		if (geosphere) geoObject = &model->boundingSphere;
		else geoObject = &model->boundingCube;

		if (splitters.size > 0) splitters.updateVaoWithLines(geoObject);
		if (samples.size > 0 && toggleSampleLines) samples.updateVaoWithLines(geoObject, rst->maindir);
		if (eslLines.size > 0 && toggle4lines) eslLines.updateVaoWithLines(geoObject);
	}


	void drawRay(int xpos, int ypos) {
		alphaLines = 0.5f;
		Ray r = cam.pixRayDirection(getDrawRayPos(xpos, ypos));
		std::vector<Ray> rays = { r };
		leaf = rst->descend(r);
		clickRay.updateVaoWithLines(rays, geoObject, rst->maindir);
		viewline = true;
		updateLines();
	}


	void picking(int xpos, int ypos) {
		Ray r = cam.pixRayDirection(getDrawRayPos(xpos, ypos));
		int primindex = -1;
		float t = 0;
		selectedPrim = -1;
		model->getIntersectionEmbree(r, primindex, t);
		if (primindex >= 0) {
			std::cout << primindex << std::endl;
			Ray extremalLine;
			std::vector<glm::vec3> edges;
			std::vector<Ray> eslEdges;

			if (rst->check1Prim(primindex, extremalLine, leaf, true, 0, true, edges, eslEdges))
				eslLines.updateVaoWithLines(std::vector<Ray>{extremalLine}, geoObject, rst->maindir);

			edgeRays.makeVaoVbo(edges);
			eslSilhEdges.updateVaoWithLines(eslEdges, geoObject);
			edgelines = true;

			selectedPrim = primindex;
			changeColoring = true;
		}
	}

	void render() {
		while (!glfwWindowShouldClose(window)) {
			glfwPollEvents();
			glUseProgram(0);
			checkColorChange();
			updateCamera(cam, width, height);
			drawElements();
			// Present result to the screen
			glfwSwapBuffers(window);
		}
	}

	// Mouse button handle function
	void mouseButtonHandler(GLFWwindow* window, int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_2 && action == GLFW_RELEASE) {
			double xpos, ypos;
			//getting cursor position
			glfwGetCursorPos(window, &xpos, &ypos);
			drawRay(xpos, ypos);

		}
		else if (button == GLFW_MOUSE_BUTTON_1 && action == GLFW_RELEASE) {
			double xpos, ypos;
			//getting cursor position
			glfwGetCursorPos(window, &xpos, &ypos);
			if (togglePicking) {
				picking(xpos, ypos);
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

	// Key handle function
	void keyboardHandler(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		//cameraKeyboardHandler(key, action);
		if (action == GLFW_PRESS) {

			switch (key)
			{
			case GLFW_KEY_1:
				leafChange(1);
				break;
			case GLFW_KEY_0:
				leafChange(-1);
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
				updateLines();
				break;
			case GLFW_KEY_6:
				toggle4lines = !toggle4lines;
				updateLines();
				//update4lines = true;
				break;
			case GLFW_KEY_7:
				togglePicking = !togglePicking;
				break;
			case GLFW_KEY_8:
				toggleGeo = !toggleGeo;
				changeGeoObject(toggleGeo);
				break;
			case GLFW_KEY_9:
				toggleLines = !toggleLines;
				break;
			default:
				break;
			}
		}
	}

	// OpenGL debug callback
	void APIENTRY debugCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam) {
		if (severity != GL_DEBUG_SEVERITY_NOTIFICATION) {
			std::cerr << "OpenGL: " << message << std::endl;
		}
	}

	void setupLights() {
		float lightdis = 3 * model->radius;
		glm::vec3 lightPos1 = lights.pos1 * lightdis;
		glm::vec3 lightPos2 = lights.pos2 * lightdis;
		glm::vec3 lightPos3 = lights.pos3 * lightdis;;

		glUseProgram(mainShader.index);
		glUniform3fv(glGetUniformLocation(mainShader.index, "lightPos1"), 1, glm::value_ptr(lightPos1));
		glUniform3fv(glGetUniformLocation(mainShader.index, "lightPos2"), 1, glm::value_ptr(lightPos2));
		glUniform3fv(glGetUniformLocation(mainShader.index, "lightPos3"), 1, glm::value_ptr(lightPos3));
	}

	void setupRenderer() {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	void setupCamera() {
		cam = PreviewCamera(model->radius * 4, model->center, 100.f, 0.1f);
		cam.aspect = width / (float)height;
	}

	void makeObjects() {
		model->boundingBox.vaoGeneration();
		model->boundingCube.vaoGeneration();
		model->boundingSphere.vaoGeneration(10, 20);

		////////////////// Objects to intersect lines with
		if (toggleGeo) geoObject = &model->boundingSphere;
		else geoObject = &model->boundingCube;

		splitters.setupVao();
		samples = LineModel(true);
		samples.setupVao();
		clickRay.setupVao(); 
		edgeRays.setupVao();
		eslSilhEdges.setupVao();
		eslLines.setupVao();

		splitters.updateVaoWithLines(rst->splitters, geoObject, rst->maindir);
	}

	void setToggles() {
		togglePoints = false;
		toggle4lines = false;
		toggleSphere = false;
		toggleBbox = false;
		toggleLines = true;
		toggleSampleLines = false;
		togglePicking = false;
		toggleGeo = false;

		leafnum = -1;
		showing = false;
		viewline = false;
		edgelines = false;
		alphaLines = 1.f;
		selectedPrim = -1;

	}

	void setControls() {
		glfwSetKeyCallback(window, keyboardHandler);
		glfwSetMouseButtonCallback(window, mouseButtonHandler);
		glfwSetCursorPosCallback(window, cursorPosHandler);
		glfwSetScrollCallback(window, scrollHandler);
	}

	void makeWindow() {
		if (!glfwInit()) std::cerr << "Failed to initialize GLFW!" << std::endl;

		//////////////////// Create window and OpenGL 4.3 debug context
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
		glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

		window = glfwCreateWindow(width, height, "Ray Space Tree", nullptr, nullptr);
		if (!window) {
			std::cerr << "Failed to create OpenGL context!" << std::endl;
			std::cout << "Press enter to close."; getchar();
		}
		glfwMakeContextCurrent(window);
	}

	void createOpenGLContext() {
		 // Initialize GLEW extension loader
		 glewExperimental = GL_TRUE;
		 glewInit();

		 // Set up OpenGL debug callback
		 glDebugMessageCallback(debugCallback, nullptr);
	}

	void initialize() {
		mainShader = Shader("standard.vert", "standard.frag");
		lineShader = Shader("line.vert", "line.frag");
		setToggles();
		//makeWindow();
		setControls();
		makeObjects();
		setupCamera();
		setupLights();
		setupRenderer();
	}

	void cleanup() {
		glfwDestroyWindow(window);
		glfwTerminate();
	}

	void visualize(int w, int h, Model* m, RaySpaceTree* r, VisComponents& vis, GLFWwindow* glfww) {// : width(w), height(h), model(m), rst(rst), visComp(vis) {
		visComp = vis;
		width = w;
		height = h;
		model = m;
		window = glfww;
		rst = r;
		initialize();
		render();
		cleanup();
	}
};