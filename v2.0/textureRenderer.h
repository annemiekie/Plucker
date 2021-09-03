#ifndef TEXTURERENDERER_H
#define TEXTURERENDERER_H

// Library for OpenGL function loading
// Must be included before GLFW
#define GLEW_STATIC
#include <GL/glew.h>
#include "shader.h"
#include "camera.h"
#include "orthocamera.h"
#include "perscamera.h"
#include "model.h"

class TextureRenderer {
public:
	TextureRenderer() {};
	Shader shader;
	GLuint rsttex;
	GLuint depthtex;
	GLuint framebuffer;
	int width;
	int height;
	Model* model;
	glm::ivec2 res;

	TextureRenderer(Shader rstshader, int rstwidth, int rstheight, Model* model) : shader(rstshader), width(rstwidth), height(rstheight), model(model) {	
		setupTexture();
		res = glm::ivec2(width, height);
	};

	void setupTexture() {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
		rsttex = generateTexture();
		depthtex = generateDepthTexture();
		framebuffer = generateFrameBuffer();
		//GLuint depthbuffer = Buffers::generateDepthBuffer(rstwidth, rstheight);
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) return;

		// Bind the shader
		glUseProgram(shader.index);

		// Bind vertex data & buffer
		glBindVertexArray(model->vao);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

		// Set viewport size
		glViewport(0, 0, width, height);
		glEnable(GL_DEPTH_TEST);
		//glDepthFunc(GL_LESS);
		glEnable(GL_CULL_FACE);
		//glCullFace(GL_BACK);
	}

	GLfloat* render(Camera * cam) {
		glm::mat4 mvp = cam->vpMatrix();
		glUniformMatrix4fv(glGetUniformLocation(shader.index, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

		glActiveTexture(GL_TEXTURE0);
		//glBindTexture(GL_TEXTURE_2D, secondtex);
		glBindTexture(GL_TEXTURE_2D, rsttex);

		glClearDepth(1.0f);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Execute draw command
		glDrawArrays(GL_TRIANGLES, 0, model->vertices.size());


		int size = width*height;
		GLfloat* pixels = new GLfloat[size];
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, pixels);

		return pixels;
	}


private:
	GLuint generateDepthBuffer() {
		GLuint depthrenderbuffer;
		glGenRenderbuffers(1, &depthrenderbuffer);
		glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
		glBindRenderbuffer(GL_RENDERBUFFER, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);
		return depthrenderbuffer;
	}

	GLuint generateFrameBuffer() {
		GLuint framebuffer;
		glGenFramebuffers(1, &framebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

		// Set "renderedTexture" as our colour attachement #0
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rsttex, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,depthtex, 0);

		GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
		glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers

		return framebuffer;
	}

	GLuint generateTexture() {
		GLuint texRST;

		glGenTextures(1, &texRST);
		glBindTexture(GL_TEXTURE_2D, texRST);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);

		//// Set interpolation for texture sampling (GL_NEAREST for no interpolation)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);

		return texRST;
	}


	GLuint generateDepthTexture() {
		GLuint texDepth;

		glGenTextures(1, &texDepth);
		glBindTexture(GL_TEXTURE_2D, texDepth);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

		//// Set interpolation for texture sampling (GL_NEAREST for no interpolation)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);

		return texDepth;
	}
};
#endif