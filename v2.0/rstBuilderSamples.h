#pragma once
#include <queue>
#include "nodeSamples.h"
#include "orthocamera.h"
#include "textureRenderer.h"
#include "rstBuilder.h"
#include "sphereSampler.h"
#include "shader.h"
#include "camera.h"
#include "rst.h"

class RSTBuilderSamples : public RSTBuilder<RSTBuilderSamples> {
private:
	//std::queue<nodeSamples> toProcess = std::queue<nodeSamples>();
	//std::vector<Ray> raytraceRays;

public:

	static void build(BuildOptions& options, RaySpaceTree* rst, VisComponents& visComp) {
		// spheresampler for now, will be cubesampling
		SphereSampler sampler = SphereSampler(&(rst->model->boundingSphere), options.noSamples);
		if (rst->alldir) sampler.createSamplesFull();
		else sampler.createSamples(rst->maindir);
		visComp.sampling_vao = sampler.vaoGeneration();
		visComp.sampling_size = sampler.samples.size();

		Orthocamera ocam(rst->model->radius, rst->model->radius, rst->model->radius * 2 > 30.f ? rst->model->radius * 2 : 30.f);
		Camera* cam = &ocam;

		// Only necessary for rasterization
		if (options.rasterizationSampling) {
			Shader rstProgram = Shader("rst.vert", "rst.frag");
			TextureRenderer texrender = TextureRenderer(rstProgram, options.width, options.height, rst->model);
			fillRST(rst, cam, texrender, sampler, options.construct, options.storeRays);
		}
		else {
			fillRSTembree(rst, cam, sampler, options);
		}

	}

	static void makeSample(glm::ivec2 res, Camera* cam, GLfloat* pixels, RaySpaceTree* rst, bool storeRays) {

		for (int y = 0; y < res.y; y++) {
			int yind = (res.y - 1 - y);
			for (int x = 0; x < res.x; x++) {
				const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
				Ray ray = cam->pixRayDirection(pixpos);
				float tri = pixels[yind * res.x + x];
				if (tri > 0) {
					float t;
					int tri_id = int(tri - 1);
					if (rst->model->getIntersectionWithPrim(tri_id, ray, t))
						rst->putPrimitive(ray, tri_id, storeRays);
					else if (rst->model->getIntersectionEmbree(ray, tri_id, t))
						rst->putPrimitive(ray, tri_id, storeRays);

				}
			}
		}
	}

	static void makeSample(glm::ivec2 res, Camera* cam, GLfloat* pixels, RaySpaceTree* rst, std::vector<std::pair<int, int>>& samples,
		std::set<int>& tris, int& count) {

		for (int y = 0; y < res.y; y++) {
			int yind = (res.y - 1 - y);
			for (int x = 0; x < res.x; x++) {
				const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
				Ray ray = cam->pixRayDirection(pixpos);
				float tri = pixels[yind * res.x + x];
				int tri_id = int(tri - 1);
				if (tri_id >= 0) {
					float t;
					if (rst->model->getIntersectionWithPrim(tri_id, ray, t))
						tris.insert(tri_id);
					else if (rst->model->getIntersectionEmbree(ray, tri_id, t))
						if (tri_id >= 0) tris.insert(tri_id);
					samples.push_back({ count, tri_id });
				}
				else if (tri_id == -1) samples.push_back({ count, -1 });
				count++;
			}
		}
	}

	static void fillMoreSamples(RaySpaceTree* rst, Camera* cam, TextureRenderer& texrender, SphereSampler& sampler) {

		std::cout << "Rasterizing camera samples to fill tree..." << std::endl;
		for (int i = 0; i < sampler.samples.size(); i++) {
			if (i % sampler.ratio == 0) continue;
			cam->setPositionAndForward(sampler.samples[i], rst->model->center);
			GLfloat* pixels = texrender.render(cam);
			makeSample(texrender.res, cam, pixels, rst, false);
		}
		int noSamples = (sampler.samples.size() - (sampler.samples.size() / sampler.ratio + 1)) * texrender.height * texrender.width;
		//std::cout << "Put " << noSamples << " samples in RST in " << diff << " ms." << std::endl;
	}

	static void fillRST(RaySpaceTree* rst,  Camera* cam, TextureRenderer& texrender, SphereSampler& sampler, int option, bool storeRays) {

		std::cout << "Rasterizing camera samples to create tree..." << std::endl;

		for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
			cam->setPositionAndForward(sampler.samples[i], rst->model->center);
			GLfloat* pixels = texrender.render(cam);
			makeSample(texrender.res, cam, pixels, rst, storeRays);
		}
	}

	static void fillRSTembree(RaySpaceTree* rst, Camera* cam, SphereSampler& sampler, BuildOptions& options) {

		std::cout << "Raytracing camera samples to create tree..." << std::endl;

		for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
			cam->setPositionAndForward(sampler.samples[i], rst->model->center);
			for (int y = 0; y < options.height; y++) {
				int yind = (options.height - 1 - y);
				for (int x = 0; x < options.width; x++) {
					const glm::vec2 pixpos{ (x + .5f) / options.width * 2.0f - 1.0f, 1.0f - (y + .5f) / options.height * 2.0f };
					Ray ray = cam->pixRayDirection(pixpos);
					if (!rst->model->boundingCube.intersectSide(rst->maindir, ray)) continue;
					int tri_id = -1;
					float t = 0.f;
					if (rst->model->getIntersectionEmbree(ray, tri_id, t, true))
						rst->putPrimitive(ray, tri_id, options.storeRays);
					//else rst->putPrimitive(ray, tri_id, storeRays, false);
				}
			}
		}
	}

	static void makeAdaptiveRST(RaySpaceTree* rst, Camera* cam, TextureRenderer& texrender, SphereSampler& sampler, int option, bool print) {
		int imgsize = texrender.width * texrender.height;
		std::vector<std::pair<int, int>> samples = std::vector<std::pair<int, int>>();
		std::set<int> tris = std::set<int>();
		std::vector<Camera*> cams = std::vector<Camera*>();
		int count = 0;
		int zerocount = 0;

		std::cout << "Rasterizing camera samples to create tree..." << std::endl;
		std::vector<Ray> raytraceRays;
		// prepare samples
		for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
			Camera* cam1 = cam->makeCopy();
			cam1->setPositionAndForward(sampler.samples[i], rst->model->center);

			GLfloat* pixels = texrender.render(cam1);
			makeSample(texrender.res, cam1, pixels, rst, samples, tris, count);
			cams.push_back(cam);
		}
		//std::cout << "Made and stored " << samples.size() << " samples in " << diff << " ms." << std::endl;
		//std::cout << "Removed " << texrender.width * texrender.height * (sampler.samples.size() / sampler.ratio + 1) - samples.size() << " samples for being too close to a primitive edge." << std::endl;

		std::cout << "Constructing the tree..." << std::endl;
		if (option == 0) rst->constructAdaptive(texrender.res, cams, samples, tris, print);
		else if (option == 1) rst->constructSmartRandom(texrender.res, cams, samples, tris);
	}

	//void compareIntersectionMethods(Camera* cam, Model* model, TextureRenderer& texrender) {
	//	/////////////// compare with raytrace
	//	//GLfloat* pixels = texrender.render(cam);

	//	for (int y = 0; y < texrender.height; y++) {
	//		int yind = (texrender.height - 1 - y);
	//		for (int x = 0; x < texrender.width; x++) {
	//			const glm::vec2 pixpos{ (x + .5f) / texrender.width * 2.0f - 1.0f, 1.0f - (y + .5f) / texrender.width * 2.0f };
	//			Ray ray = cam->pixRayDirection(pixpos);
	//			//int triRast = int(pixels[yind * texrender.width + x])-1;
	//			float t1 = 0.f; float t2 = 0.f;
	//			int triTracePlucker = -1;
	//			model->getIntersectionNoAcceleration(ray, triTracePlucker, t1);
	//			int triTraceEmbree = -1;
	//			model->getIntersectionEmbree(ray, triTraceEmbree, t2);
	//			//if (triRast != triTracePlucker || triRast != triTraceEmbree || 
	//			if (triTraceEmbree != triTracePlucker) {
	//				//std::cout << "Rasterization: " << triRast << 
	//				std::cout << "Trace Plucker: " << triTracePlucker << " Trace Embree: " << triTraceEmbree << std::endl;
	//			}
	//			if (triTraceEmbree >= 0 && std::fabs(t1 - t2) > 1E-8) {
	//				std::cout << "Trace Plucker depth: " << t1 << "Trace Embree depth: " << t2 << std::endl;
	//			}

	//		}
	//	}
	//}

};

