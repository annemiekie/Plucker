#pragma once
// Library for OpenGL function loading
// Must be included before GLFW
#define GLEW_STATIC
#include <GL/glew.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

class Shader {
public:
	GLuint index;


	Shader(std::string vert, std::string frag) {
		////////////////// Load and compile rst shader program
		index = glCreateProgram();

		std::string vertexShaderCode = readFile(vert);
		const char* vertexShaderCodePtr = vertexShaderCode.data();

		GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertexShader, 1, &vertexShaderCodePtr, nullptr);
		glCompileShader(vertexShader);

		std::string fragmentShaderCode = readFile(frag);
		const char* fragmentShaderCodePtr = fragmentShaderCode.data();

		GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragmentShader, 1, &fragmentShaderCodePtr, nullptr);
		glCompileShader(fragmentShader);

		if (!checkShaderErrors(vertexShader) || !checkShaderErrors(fragmentShader)) {
			std::cerr << "Shader(s) failed to compile!" << std::endl;
			return;
		}

		// Combine vertex and fragment shaders into single shader program
		glAttachShader(index, vertexShader);
		glAttachShader(index, fragmentShader);
		glLinkProgram(index);

		if (!checkProgramErrors(index)) {
			std::cerr << "Shadow program failed to link!" << std::endl;
			return;
		}
	};

private:
	// Helper function to read a file like a shader
	std::string readFile(const std::string& path) {
		std::ifstream file(path, std::ios::binary);

		std::stringstream buffer;
		buffer << file.rdbuf();

		return buffer.str();
	};

	bool checkShaderErrors(GLuint shader) {
		// Check if the shader compiled successfully
		GLint compileSuccessful;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compileSuccessful);

		// If it didn't, then read and print the compile log
		if (!compileSuccessful) {
			GLint logLength;
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);

			std::vector<GLchar> logBuffer(logLength);
			glGetShaderInfoLog(shader, logLength, nullptr, logBuffer.data());

			std::cerr << logBuffer.data() << std::endl;

			return false;
		}
		else {
			return true;
		}
	};


	bool checkProgramErrors(GLuint program) {
		// Check if the program linked successfully
		GLint linkSuccessful;
		glGetProgramiv(program, GL_LINK_STATUS, &linkSuccessful);

		// If it didn't, then read and print the link log
		if (!linkSuccessful) {
			GLint logLength;
			glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);

			std::vector<GLchar> logBuffer(logLength);
			glGetProgramInfoLog(program, logLength, nullptr, logBuffer.data());

			std::cerr << logBuffer.data() << std::endl;

			return false;
		}
		else {
			return true;
		}
	};




};