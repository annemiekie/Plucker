#ifndef SAMPLER_H
#define SAMPLER_H

#include <GL/glew.h>
// Library for vertex and matrix math
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
//#include <vector>
#include "sphere.h"
#include "cube.h"
#include "sampler.h"
#include "orthocamera.h"

struct SphereSampler : Sampler {

	Sphere* sphere;
	int camLoc = 0;
	Orthocamera currentCam;
	glm::vec3 camLookat;

	SphereSampler() : Sampler() {};

	SphereSampler(Sphere* sphere, int N, int wh, glm::vec3 lookat, glm::vec3 maindir, bool alldir, int ratio = 1) : 
		sphere(sphere), camLookat(lookat), Sampler(N, wh, ratio, maindir, alldir) {
		currentCam = Orthocamera(sphere->radius, sphere->radius, sphere->radius * 2 > 30.f ? sphere->radius * 2 : 30.f);
		createMainSamples();
		currentCam.setPositionAndForward(main_grid[0], lookat);
	};

	static glm::vec2 getPixPosFromIndex(int w, int h, int x, int y) {
		return { (x + .5f) / w * 2.0f - 1.0f, 1.0f - (y + .5f) / h * 2.0f };
	}

	void getIndices(int raynr, int& camNr, int& x, int& y) {
		camNr = raynr / detailGridNum;
		int rayInCam = raynr % detailGridNum;
		y = detailGridDim - rayInCam / detailGridDim;
		x = rayInCam % detailGridDim;
	}

	virtual Ray getSample(int raynr) override {
		int camNr, y, x;
		getIndices(raynr, camNr, x, y);
		const glm::vec2 pixpos = getPixPosFromIndex(detailGridDim, detailGridDim, x, y);
		Orthocamera cam(detailGridDim, detailGridDim);
		cam.setPositionAndForward(main_grid[camNr], camLookat);
		return cam.pixRayDirection(pixpos);
	}

	virtual Ray getNextSample(bool inverseRatio) override {
		if (inverseRatio && counter % ratio == 0) counter++;
		if (counter % (detailGridDim * detailGridDim) == 0)
			currentCam.setPositionAndForward(main_grid[counter / detailGridNum], camLookat);
		int camNr, y, x;
		getIndices(counter, camNr, x, y);
		counter += inverseRatio ? 1 : ratio;
		return currentCam.pixRayDirection(getPixPosFromIndex(detailGridDim, detailGridDim, x, y));
	}

	virtual void createSamplesFull() override {
		double x, y, z;
		double factorTheta = 2. * glm::pi<double>() * 2. / (1. + sqrt(5));

		for (int i = 0; i <= mainGridDim; i++) {
			double cosphi = 1. - 2. * i / mainGridDim;//sgn*
			double theta = i * factorTheta;//sgn*
			double phi = acos(cosphi);
			x = sphere->radius * cos(theta) * sin(phi);
			y = sphere->radius * sin(theta) * sin(phi);
			z = sphere->radius * cos(phi);

			glm::dvec3 sample;
			sample = glm::dvec3(x, y, z);
			main_grid.push_back(sample + sphere->center);
		}
	}

	void findMinMax(glm::vec3& min, glm::vec3& max) {
		min = { 1000, 1000, 1000 };
		max = { -1000, -1000, -1000 };
		for (glm::dvec3& s : main_grid) {
			if (s.x < min.x) min.x = s.x;
			else if (s.x > max.x) max.x = s.x;
			if (s.y < min.y) min.y = s.y;
			else if (s.y > max.y) max.y = s.y;
			if (s.z < min.z) min.z = s.z;
			else if (s.z > max.z) max.z = s.z;
		}

	}

	virtual void createSamples() override {
		int sgn = glm::dot(glm::vec3(-1), maindir);
		float phi, theta, sintheta;
		float x, y, z;
		float factorPhi = 2.f * glm::pi<float>() * 2.f / (1.f + sqrtf(5));
		float factorTheta = 2.f / (2.f * mainGridDim + 1.f);
		glm::vec3 dir;
		if (fabsf(maindir.x)) dir = glm::vec3(1, 2, 0);
		else if (fabsf(maindir.y)) dir = glm::vec3(2, 0, 1);
		else if (fabsf(maindir.z)) dir = glm::vec3(0, 1, 2);
		else return;

		for (int i = 0; i <= mainGridDim; i++) {
			phi = sgn * i * factorPhi;//
			sintheta = sgn * i * factorTheta;//
			theta = asin(sintheta);
			x = sphere->radius * cos(theta) * cos(phi);
			y = sphere->radius * cos(theta) * sin(phi);
			z = sphere->radius * sintheta;

			//if (fabsf(z) < fabsf(x) || fabsf(z) < fabsf(y)) continue;

			glm::dvec3 sample;
			sample[dir.x] = x;
			sample[dir.y] = y;
			sample[dir.z] = z;
			main_grid.push_back(sample + sphere->center);
		}
	}



};

#endif