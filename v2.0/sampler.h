#pragma once
#include "ray.h"

struct Sampler {
	std::vector<glm::dvec3> main_grid;
	uint64_t mainGridDim = 0;
	uint64_t detailGridDim = 0;
	uint64_t detailGridNum = 0;
	uint64_t counter = 0;
	uint64_t ratio = 1;
	uint64_t nrOfSamples = 0;
	glm::vec3 maindir;
	bool alldir;
	GLuint vao = 0;


	Sampler() {};

	Sampler(int N_main, int N_lower, int ratio, glm::vec3 maindir, bool alldir)
		: mainGridDim(N_main), detailGridDim(N_lower), ratio(ratio), maindir(maindir), alldir(alldir) {
		detailGridNum = detailGridDim * detailGridDim;
	};

	void createMainSamples() {
		if (alldir) createSamplesFull();
		else 		createSamples();
		nrOfSamples = main_grid.size() * detailGridDim * detailGridDim;
	};

	virtual Ray getNextSample(bool inverseRatio) {
		return Ray();
	};

	virtual Ray getSample(uint64_t raynr) {
		return Ray();
	};

	virtual void createSamples() {};
	virtual void createSamplesFull() {};

	int vaoGeneration() {
		GLuint pointVBO;
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glGenBuffers(1, &pointVBO);
		glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * main_grid.size(), &main_grid[0].x, GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (void*)0);
		return vao;
	}
};