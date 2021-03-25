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
	GLuint vao = 0;
	GLuint vao2 = 0;

	std::vector<Vertex> vertices = std::vector<Vertex>();
	std::vector<glm::vec3> vertices2 = std::vector<glm::vec3>();
	std::vector<int> indices = std::vector<int>();
	std::vector<std::vector<int>> triPerVertex;

	int primsize = 0;
	Cube boundingBox;
	glm::vec3 center = glm::vec3(0, 0, 0);
	float radius = 0.f;

	Model(const char* filename) {
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string err;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename)) {
			std::cerr << err << std::endl;
			return;
		}

		int ind = 0;
		int count = 0;
		float n = 1000.f;
		float maxx = -n, maxy = -n, maxz = -n;
		float minx = n, miny = n, minz = n;

		vertices2 = std::vector<glm::vec3>(attrib.vertices.size() / 3);
		indices = std::vector<int>();
		triPerVertex = std::vector<std::vector<int>>(attrib.vertices.size() / 3);

		for (int i = 0; i < attrib.vertices.size(); i+=3) {
			glm::vec3 vertex = {
				attrib.vertices[i],
				attrib.vertices[i + 1],
				attrib.vertices[i + 2]
			};
			vertices2[i/3] = vertex;
		}

		std::vector<glm::vec3> centroid(3);
		for (const auto& shape : shapes) {
			primsize += shape.mesh.indices.size() / 3;

			for (const auto& index : shape.mesh.indices) {
				indices.push_back(index.vertex_index);
				Vertex vertex = {};

				// Retrieve coordinates for vertex by index
				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2]
				};

				vertex.bary = { count % 3 == 0, (count + 1) % 3 == 0, (count + 2) % 3 == 0 };


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

				triPerVertex[index.vertex_index].push_back(vertex.id - 1);

				vertices.push_back(vertex);
			}
		}

		for (int i = 0; i < vertices.size(); i+=3) {
			glm::vec3 centroid = (vertices[i].pos + vertices[i+1].pos + vertices[i+2].pos) / 3.f;
			vertices[i].center = centroid;
			vertices[i + 1].center = centroid;
			vertices[i + 2].center = centroid;
		}

		boundingBox = Cube(glm::vec3(minx, miny, minz), glm::vec3(maxx, maxy, maxz));
		center = glm::vec3((minx + maxx) / 2.f, (miny + maxy) / 2.f, (minz + maxz) / 2.f);
		radius = glm::length(glm::vec3(maxx, maxy, maxz) - center);

		createVAO();
		createIBO();
	};
	
	void createIBO() {
		glGenVertexArrays(1, &vao2); // make VAO
		glBindVertexArray(vao2);

		// copy interleaved vertex data (V/N/T) to VBO
	   // GLuint vbo[2];
		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, vertices2.size() * sizeof(glm::vec3), &vertices2[0], GL_STATIC_DRAW);

		// copy index data to VBO
		GLuint ibo;
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);

		// The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
		// Stride is the distance in bytes between vertices
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
		glEnableVertexAttribArray(0);

	}

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

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, normal)));
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, bary)));
		glEnableVertexAttribArray(2);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, center)));
		glEnableVertexAttribArray(3);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, selected)));
		glEnableVertexAttribArray(4);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, id)));
		glEnableVertexAttribArray(5);
	}

};

#endif