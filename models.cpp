#include "models.h"


GLuint Models::vaoLineGeneration(std::vector<glm::vec3>& lines) {

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


GLuint Models::quadGeneration() {
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



GLuint Models::vaoLineGenerationWithColor(std::vector<glm::vec3>& lines, std::vector<glm::vec3>& colorpos, Cube* cube) {

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