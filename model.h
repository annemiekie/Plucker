#ifndef MODEL_H
#define MODEL_H

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

class Model {
public:
	GLuint vao;
	std::vector<Vertex>& vertices = std::vector<Vertex>();
	int primsize = 0;
	Cube boundingBox;
	glm::vec3 center = glm::vec3(0, 0, 0);
	float radius = 0.f;

	Model(const char* filename, int shape) {
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string err;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename)) {
			std::cerr << err << std::endl;
			return;
		}

		primsize = shapes[shape].mesh.indices.size() / 3;
		int ind = 0;
		int count = 0;
		float n = 1000.f;
		float maxx = -n, maxy = -n, maxz = -n;
		float minx = n, miny = n, minz = n;
		// Read triangle vertices from OBJ file
		//for (const auto& shape : shapes) {
		for (const auto& index : shapes[shape].mesh.indices) {
			Vertex vertex = {};

			// Retrieve coordinates for vertex by index
			vertex.pos = {
				attrib.vertices[3 * index.vertex_index + 0],
				attrib.vertices[3 * index.vertex_index + 1],
				attrib.vertices[3 * index.vertex_index + 2]
			};

			if (vertex.pos.x > maxx) maxx = vertex.pos.x;
			else if (vertex.pos.x < minx) minx = vertex.pos.x;
			if (vertex.pos.y > maxy) maxy = vertex.pos.y;
			else if (vertex.pos.y < miny) miny = vertex.pos.y;
			if (vertex.pos.z > maxz) maxz = vertex.pos.z;
			else if (vertex.pos.z < minz) minz = vertex.pos.z;

			// Retrieve components of normal by index
			vertex.normal = {
				attrib.normals[3 * index.normal_index + 0],
				attrib.normals[3 * index.normal_index + 1],
				attrib.normals[3 * index.normal_index + 2]
			};

			if (count % 3 == 0) {
				ind++;
			}
			vertex.id = (1.f * ind);
			count++;

			vertices.push_back(vertex);
		}
		boundingBox = Cube(glm::vec3(minx, miny, minz), glm::vec3(maxx, maxy, maxz));
		center = glm::vec3((minx + maxx) / 2.f, (miny + maxy) / 2.f, (minz + maxz) / 2.f);
		radius = glm::length(glm::vec3(maxx, maxy, maxz) - center);

		createVAO();
	};


	void createVAO() {
		//////////////////// Create Vertex Buffer Object
		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

		// Bind vertex data to shader inputs using their index (location)
		// These bindings are stored in the Vertex Array Object
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
	}

};

#endif