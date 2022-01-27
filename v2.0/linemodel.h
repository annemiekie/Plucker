#pragma once

#ifndef LINEMODEL_H
#define LINEMODEL_H

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

class LineModel {
public:
	GLuint vao;
	GLuint vbo;
	GLuint vboc;
	bool color = false;
	int size = 0;
	std::vector<Ray> rays;

	LineModel(bool color = false) : color(color) {};

	void setupVao() {
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
		glEnableVertexAttribArray(0);
		if (color) {
			glGenBuffers(1, &vboc);
			glBindBuffer(GL_ARRAY_BUFFER, vboc);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
			glEnableVertexAttribArray(1);
		}
	}

	void updateVaoWithLines(std::vector<Ray> rays_new, GeoObject* object, glm::vec3 maindir = glm::vec3(0)) {
		if (rays_new.size() == 0) return;
		rays = rays_new;
		updateVaoWithLines(object, maindir);
	}

	void updateVaoWithLines(GeoObject* object, glm::vec3 maindir = glm::vec3(0)) {
		std::vector<glm::vec3> lines, colors;
		getLinesInGeo(object, lines, colors, maindir);
		if (lines.size()) makeVaoVbo(lines, colors);
	}

	void makeVaoVbo(std::vector<glm::vec3>& lines, std::vector<glm::vec3>& colors = std::vector<glm::vec3>()) {
		size = lines.size();
		glBindVertexArray(vao);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * lines.size(), &lines[0].x, GL_STATIC_DRAW);

		if (color) {
			glBindBuffer(GL_ARRAY_BUFFER, vboc);
			glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * colors.size(), &colors[0].x, GL_STATIC_DRAW);
		}
	}


	void getLinesInGeo(GeoObject* object, std::vector<glm::vec3>& lines, std::vector<glm::vec3>& colors = std::vector<glm::vec3>(), glm::vec3 maindir = glm::vec3(0)) {
		glm::vec3 line, col;
		for (Ray& r : rays) {
			if (r.direction.x == 0 && r.direction.y == 0 && r.direction.z == 0) r.get3DfromPlucker();
			float s, t;
			glm::vec3 c;
			if (object->intersect(r, s, t, color, maindir, c)) {
				lines.push_back(r.origin + (double)s * r.direction);
				lines.push_back(r.origin + (double)t * r.direction);
				if (color) {
					colors.push_back(c);
					colors.push_back(c);
				}
			}
		}
	}

	void draw() {
		glBindVertexArray(vao);
		glDrawArrays(GL_LINES, 0, size);
	}

};

#endif

