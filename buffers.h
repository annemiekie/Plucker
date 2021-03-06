#ifndef BUFFERS_H
#define BUFFERS_H

// Library for OpenGL function loading
// Must be included before GLFW
#define GLEW_STATIC
#include <GL/glew.h>


namespace Buffers {


	GLuint generateDepthBuffer(const int rstwidth, const int rstheight) {
		GLuint depthrenderbuffer;
		glGenRenderbuffers(1, &depthrenderbuffer);
		glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, rstwidth, rstheight);
		glBindRenderbuffer(GL_RENDERBUFFER, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);
		return depthrenderbuffer;
	}

	GLuint generateFrameBuffer(GLuint texRST, GLuint texDepth) {
		GLuint framebuffer;
		glGenFramebuffers(1, &framebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

		// Set "renderedTexture" as our colour attachement #0
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texRST, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, texDepth, 0);

		GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
		glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers

		return framebuffer;
	}

	GLuint generateTexture(const int rstwidth, const int  rstheight) {
		GLuint texRST;

		glGenTextures(1, &texRST);
		glBindTexture(GL_TEXTURE_2D, texRST);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, rstwidth, rstheight, 0, GL_RED, GL_FLOAT, nullptr);

		//// Set interpolation for texture sampling (GL_NEAREST for no interpolation)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);

		return texRST;
	}


	GLuint generateDepthTexture(const int rstwidth, const int  rstheight) {
		GLuint texDepth;

		glGenTextures(1, &texDepth);
		glBindTexture(GL_TEXTURE_2D, texDepth);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, rstwidth, rstheight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

		//// Set interpolation for texture sampling (GL_NEAREST for no interpolation)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);

		return texDepth;
	}
}
#endif