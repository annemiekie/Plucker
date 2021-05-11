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
#include "geoObject.h"
#include "ray.h"

namespace Models {

	GLuint quadGeneration();

};


#endif