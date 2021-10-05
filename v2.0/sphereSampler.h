#ifndef SAMPLER_H
#define SAMPLER_H

#include <GL/glew.h>
// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include "Sphere.h"

struct SphereSampler {
	std::vector<glm::vec3> samples = std::vector<glm::vec3>();
	Sphere* sphere;
	GLuint vao = 0;
	int ratio = 1;
	int N = 0;
	int camLoc = 0;
	
	SphereSampler(Sphere* sphere, int N, int ratio) : sphere(sphere), N(N), ratio(ratio) {};

	void createSamplesFull() {
		float x, y, z;
		float factorTheta = 2.f * glm::pi<float>() * 2.f / (1.f + sqrtf(5));

		for (int i = 0; i <= N; i++) {
			float cosphi = 1.f - 2.f*i/N;//sgn*
			float theta = i * factorTheta;//sgn*
			float phi = acos(cosphi);
			x = sphere->radius * cos(theta) * sin(phi);
			y = sphere->radius * sin(theta) * sin(phi);
			z = sphere->radius * cos(phi);

			glm::vec3 sample;
			sample = glm::vec3(x, y, z);
			samples.push_back(sample + sphere->center);
		}
	}

	void findMinMax(glm::vec3 &min, glm::vec3 &max) {
		min = { 1000, 1000, 1000 };
		max = { -1000, -1000, -1000 };
		for (glm::vec3& s : samples) {
			if (s.x < min.x) min.x = s.x;
			else if (s.x > max.x) max.x = s.x;
			if (s.y < min.y) min.y = s.y;
			else if (s.y > max.y) max.y = s.y;
			if (s.z < min.z) min.z = s.z;
			else if (s.z > max.z) max.z = s.z;
		}

	}

	void createSamples(int sgn, char maindir) {
		sgn *= -1;
		float phi, theta, sintheta;
		float x, y, z;
		float factorPhi = 2.f * glm::pi<float>() * 2.f / (1.f + sqrtf(5));
		float factorTheta = 2.f / (2.f * N + 1.f);
		glm::vec3 dir;
		if (maindir == 'X') dir = glm::vec3(1, 2, 0);
		else if (maindir == 'Y') dir = glm::vec3(2, 0, 1);
		else if (maindir == 'Z') dir = glm::vec3(0, 1, 2);
		else return;

		for (int i = 0; i <= N; i++) {
			phi = sgn * i * factorPhi;//
			sintheta = sgn * i * factorTheta;//
			theta = asin(sintheta);
			x = sphere->radius * cos(theta) * cos(phi);
			y = sphere->radius * cos(theta) * sin(phi);
			z = sphere->radius * sintheta;

			if (fabsf(z) < fabsf(x) || fabsf(z) < fabsf(y)) continue;

			glm::vec3 sample;
			sample[dir.x] = x;
			sample[dir.y] = y;
			sample[dir.z] = z;

			samples.push_back(sample + sphere->center);
		}
	}

	void vaoGeneration() {
		GLuint pointVBO;
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glGenBuffers(1, &pointVBO);
		glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * samples.size(), &samples[0].x, GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (void*)0);
	}


};

#endif