#ifndef MODELS_H
#define MODELS_H

#define GLEW_STATIC
#include <GL/glew.h>


#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>

//// Library for loading .OBJ model
//#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "vertex.h"
#include "cube.h"

namespace Models {

	GLuint sphere(int lats, int longs, int &num, float rad, GLuint& m_vboIndex);
	GLuint vaoLineGeneration(std::vector<glm::vec3>& lines);
	GLuint loadModel(std::vector<Vertex>& vertices, int& size, const char* filename);
	GLuint quadGeneration();
	GLuint cubelineGeneration(glm::vec3& min, glm::vec3& max);
	GLuint vaoLineGenerationWithColor(std::vector<glm::vec3>& lines, std::vector<glm::vec3>& colorpos, Cube* cube);
};


#endif