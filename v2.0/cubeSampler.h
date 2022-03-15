#include "cube.h"
#include "sampler.h"
#include "square.h"

struct CubeSampler : Sampler {

	Cube* cube;
	glm::dvec3 currentDome;
	glm::ivec3 dir;
	double sgn = 1;

	CubeSampler() : Sampler() {};


	CubeSampler(Cube* cube, int N, int wh, glm::vec3 maindir, bool alldir, int ratio = 1) :
				cube(cube), Sampler(round(sqrt(1.f*N)), wh, ratio, maindir, alldir) {

		if (fabsf(maindir.x)) dir = glm::ivec3(1, 2, 0);
		else if (fabsf(maindir.y)) dir = glm::ivec3(2, 0, 1);
		else if (fabsf(maindir.z)) dir = glm::ivec3(0, 1, 2);
		sgn = glm::dot(maindir, glm::vec3(1));
		createMainSamples();
		currentDome = main_grid[0];
	};

	//virtual Ray getSample ()

	virtual Ray getNextSample(bool inverseRatio) override {
		if (inverseRatio && counter % ratio == 0) counter++;
		if (counter % detailGridNum == 0) currentDome = main_grid[counter / detailGridNum];
		// dodecahedron?
		// or subdivide into bands
		// or just theta phi for now
		int rayInDome = counter % detailGridNum;
		int theta = rayInDome / detailGridDim;
		int phi = rayInDome % detailGridDim;
		
		double deltaTh = (.5 * glm::pi<double>()) / (1. * detailGridDim);
		double deltaPh = glm::pi<double>() / (1. * detailGridDim);

		glm::dvec3 rayDir;
		rayDir[dir.x] = cos(deltaPh * phi) * sin(deltaTh * theta);
		rayDir[dir.y] = sin(deltaPh * phi) * sin(deltaTh * theta);
		rayDir[dir.z] = cos(theta * deltaTh);
		counter += inverseRatio ? 1 : ratio;
		return Ray(currentDome, currentDome + sgn * rayDir);
	};

	virtual void createSamples() override {
		Square side = cube->getCubeSideSquare(maindir);
		glm::vec3 min = side.getMin();
		glm::vec3 max = side.getMax();

		glm::vec3 sample;
		sample[dir.z] = min[dir.z];
		for (int i = 0; i <= mainGridDim; i++) {
			sample[dir.x] = min[dir.x] + i * (max[dir.x] - min[dir.x]) / (1.*mainGridDim);
			for (int j = 0; j <= mainGridDim; j++) {
				sample[dir.y] = min[dir.y] + j * (max[dir.y] - min[dir.y]) / (1. * mainGridDim);
				main_grid.push_back(sample);
			}
		}
	}
};