#ifndef MODELS_H
#define MODELS_H

#define GLEW_STATIC
#include <GL/glew.h>


#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "vertex.h"
#include "cube.h"

namespace Models {

	GLuint vaoLineGeneration(std::vector<glm::vec3>& lines);
	GLuint quadGeneration();
	//GLuint cubelineGeneration(glm::vec3& min, glm::vec3& max);
	GLuint vaoLineGenerationWithColor(std::vector<glm::vec3>& lines, std::vector<glm::vec3>& colorpos, Cube* cube);
	GLuint rayVao(Ray& r, glm::vec3 color, Cube* cube, glm::vec3 mainDir);

};


#endif