#pragma once
#include <queue>
#include <chrono>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "nodeSamples.h"
#include "orthocamera.h"
#include "textureRenderer.h"
#include "rstBuilder.h"
#include "sphereSampler.h"
#include "cubeSampler.h"
#include "shader.h"
#include "camera.h"
#include "rst.h"
#include "sample.h"
#include "sampler.h"

class RSTBuilderSamples : public RSTBuilder<RSTBuilderSamples> {
public:

	static void build(Options::BuildOptions& options, RaySpaceTree* rst, VisComponents& visComp, bool print) {
		// spheresampler for now, will be cubesampling
		Sampler* sampler;
		SphereSampler sphereSampler;
		CubeSampler cubeSampler;
		if (options.samplingtype == Options::samplingType::UNIFORM_SPHERE) {
			sphereSampler = SphereSampler(&(rst->model->boundingSphere), options.noSamples,
				options.s_w, rst->model->center, rst->maindir, rst->alldir, options.sampleStoreFillRatio);
			sampler = &sphereSampler;
		}
		else {
			cubeSampler = CubeSampler(&(rst->model->boundingCube), options.noSamples,
				options.s_w, rst->maindir, rst->alldir, options.sampleStoreFillRatio);
			sampler = &cubeSampler;
		}
		
		visComp.sampling_vao = sampler->vaoGeneration();
		visComp.sampling_size = sampler->main_grid.size();

		//TextureRenderer texrender;
		//if (options.rasterizationSampling) {
		//	Shader rstProgram = Shader("rst.vert", "rst.frag");
		//	texrender = TextureRenderer(rstProgram, options.s_w, options.s_h);
		//}

		if (options.construct == Options::ADAPTIVE) createAdaptive(rst, options, sampler);// , texrender);
		else fill(rst, options, sampler);// , texrender);

	}

	//static glm::vec2 getPixPosFromIndex(int w, int h, int x, int y) {
	//	return { (x + .5f) / w * 2.0f - 1.0f, 1.0f - (y + .5f) / h * 2.0f };
	//}

	//static void makeSample(RaySpaceTree* rst, Options::BuildOptions& options, Sampler& sampler, bool fillMore, GLfloat* pixels = 0,
	//						std::vector<SampleInd>& samples = std::vector<SampleInd>(), int camNum = 0,
	//						std::set<int>& triangles = std::set<int>()) {
	//
	//	for (int y = 0; y < options.s_h; y++) {
	//		int yind = (options.s_h - 1 - y);
	//		for (int x = 0; x < options.s_w; x++) {
	//			Ray ray = sampler.getNextSample();
	//			//Ray ray(glm::vec3(-30, 0,-5), glm::vec3(30, 0, 0));
	//			if (!rst->alldir && !rst->model->boundingCube.intersectSide(rst->maindir, ray)) {
	//				//if (options.storeAllSamples) rst->putPrimitive(ray, -1, true, false);
	//				continue;
	//			}
	//
	//			int pixloc = yind * options.s_w + x;
	//			int tri = -1;
	//			if (options.rasterizationSampling) tri = int(pixels[pixloc]) - 1 ;
	//			else {
	//				float t = 0.f;
	//				rst->model->getIntersectionEmbree(ray, tri, t, true);
	//			}
	//			if (options.construct == Options::ADAPTIVE) {
	//				samples.push_back(SampleInd({ tri, camNum * options.s_h * options.s_w + pixloc }));
	//				if (tri >= 0) triangles.insert(tri);
	//			}
	//			else if (tri >= 0) {
	//				rst->putPrimitive(ray, tri, !fillMore && options.storeSamples);
	//			}
	//			if (options.storeAllSamples && !fillMore) {
	//				rst->putPrimitive(ray, tri, true, false);
	//			}
	//		}
	//	}
	//}

	static void makeSample(RaySpaceTree* rst, Options::BuildOptions& options, Sampler* sampler,// GLfloat* pixels = 0,
							std::vector<SampleInd>& samples = std::vector<SampleInd>(),
							std::set<int>& triangles = std::set<int>()) {

		Ray ray = sampler->getNextSample(options.fillMoreSamples);
		if (!rst->alldir && !rst->model->boundingCube.intersectSide(rst->maindir, ray)) return;
		//if (options.storeAllSamples) rst->putPrimitive(ray, -1, true, false);

		int tri = -1;
		//		if (options.rasterizationSampling) tri = int(pixels[pixloc]) - 1;
		float t = 0.f;
		rst->model->getIntersectionEmbree(ray, tri, t, true);
		if (options.construct == Options::ADAPTIVE) {
			samples.push_back(SampleInd({ tri, sampler->counter }));
			if (tri >= 0) triangles.insert(tri);
		}
		else {
			if (tri >= 0)
				rst->putPrimitive(ray, tri, !options.fillMoreSamples && options.storeSamples, true, options.exact, options.exactStartLevel);
			else if (options.storeAllSamples && !options.fillMoreSamples)
				rst->putPrimitive(ray, tri, true, false, options.exact, options.exactStartLevel);
		}


	};

	static void fill(RaySpaceTree* rst, Options::BuildOptions& options, Sampler* sampler, TextureRenderer& texrender = TextureRenderer()) {

		std::cout << "Creating camera samples to fill tree..." << std::endl;
		while (sampler->counter < sampler->nrOfSamples) {
			//GLfloat* pixels = 0;
			//if (options.rasterizationSampling) pixels = texrender.render((SphereSampler)(sampler.currentCam), rst->model);
			makeSample(rst, options, sampler);
		}
		//for (int i = 0; i < sampler.main_grid.size(); i++) {
		//	//cam->setPositionAndForward(sampler.main_grid[i], rst->model->center);
		//	//GLfloat* pixels = 0;
		//	//if (options.rasterizationSampling) pixels = texrender.render(cam, rst->model);
		//
		//	if (options.fillMoreSamples && i % sampler.ratio != 0)
		//		makeSample(rst, options, sampler, true, pixels);
		//	else 
		//		makeSample(rst, options, sampler, false, pixels);
		//}
	}

	static void createAdaptive(RaySpaceTree* rst, Options::BuildOptions& options, Sampler* sampler, TextureRenderer& texrender = TextureRenderer()) {
		
		std::vector<SampleInd>& samples = std::vector<SampleInd>();
		std::set<int> triangles = std::set<int>();
		//std::vector<Camera*> cams = std::vector<Camera*>();

		std::cout << "Rasterizing camera samples to create tree..." << std::endl;
		std::vector<Ray> raytraceRays;
		// prepare samples
		//for (int i = 0; i < sampler.main_grid.size(); i += sampler.ratio) {
			//Camera* cam1 = cam->makeCopy();
			//cam1->setPositionAndForward(sampler.main_grid[i], rst->model->center);
		//	GLfloat* pixels = 0;
		//	if (options.rasterizationSampling) pixels = texrender.render(cam1, rst->model);
		while (sampler->counter < sampler->nrOfSamples) makeSample(rst, options, sampler, samples, triangles);
		std::cout << "Constructing the tree..." << std::endl;
		constructAdaptive(rst, options, sampler, samples, triangles);
	}

	//static Ray getRayFromNum(std::vector<Camera*>& cams, int raynr,int w, int h) {
	//	int camNr = raynr / (w * h);
	//	int rayInCam = raynr % (w * h);
	//	int y = h - rayInCam / w;
	//	int x = rayInCam % w;
	//	const glm::vec2 pixpos = getPixPosFromIndex(w, h, x, y);
	//	return cams[camNr]->pixRayDirection(pixpos);
	//}

	static void constructAdaptive(RaySpaceTree* rst, Options::BuildOptions& options, Sampler* sampler, std::vector<SampleInd>& samples, std::set<int>& tris) {
		std::set<Edge*, Edge::cmp_ptr> edges;
		for (int i : tris) for (Edge* e : rst->model->triangles[i]->edges) edges.insert(e);

		std::queue<NodeSamples> toProcess;
		NodeSamples ns = { 0, samples, tris, edges, rst->rootNode };
		constructAdaptive(rst, ns, options, sampler, toProcess);
		samples = std::vector <SampleInd>();
		tris = std::set<int>();

		while (!toProcess.empty()) {
			ns = toProcess.front();
			constructAdaptive(rst, ns, options, sampler, toProcess);
			toProcess.pop();

		}
	}

	static void constructAdaptive(RaySpaceTree *rst, NodeSamples& ns, Options::BuildOptions& options, Sampler* sampler,
									 std::queue<NodeSamples>& toProcess) {
		auto start_time1 = std::chrono::high_resolution_clock::now();

		std::cout << "level: " << ns.level << std::endl;
		bool checkmore = false;
		if (ns.triangles.size() == 1) {
			for (auto s : ns.samples)
				if (s.prim < 0) {
					checkmore = true;
					break;
				}
		}

		if (ns.level <= rst->depth, (ns.triangles.size() > 1 || checkmore)) {
			std::cout << "Node with " << ns.samples.size() << " samples" << std::endl;
			std::cout << "Finding best splitter..." << std::endl;

			auto start_time = std::chrono::high_resolution_clock::now();
			Split splitter;

			//change this to all nodesamples
			NodeSamples nsLeft, nsRight;

			if (getBestSplitter(rst, splitter, options, nsLeft, nsRight, ns, sampler, 10,true)) {

				auto end_time = std::chrono::high_resolution_clock::now();
				auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
				std::cout << "Found splitter in " << diff << " ms" << std::endl;

				ns.node->splitter = splitter;
				rst->splitters.push_back(splitter);

				Node* leftnode = new Node(ns.node->index * 2 + 1, ns.level + 1);
				ns.node->leftNode = leftnode;
				leftnode->parent = ns.node;
				rst->nodes.push_back(leftnode);
				nsLeft.node = ns.node->leftNode;
				nsLeft.level = ns.level + 1;
				toProcess.push(nsLeft);
				//constructAdaptive(level + 1, node->leftNode, res, cams, samplesL, trisL, print);

				Node* rightnode = new Node(ns.node->index * 2 + 2, ns.level + 1);
				rightnode->parent = ns.node;
				ns.node->rightNode = rightnode;
				rst->nodes.push_back(rightnode);
				nsRight.node = ns.node->rightNode;
				nsRight.level = ns.level + 1;
				toProcess.push(nsRight);
				//constructAdaptive(level + 1, node->rightNode, res, cams, samplesR, trisR, print);
				return;
			}
		}
		std::cout << "Leaf with " << ns.samples.size() << " samples and " << ns.triangles.size() << " triangles" << std::endl;
		
		rst->noLeaves++;
		ns.node->leaf = true;

		if (ns.triangles.size() > 0) {
			ns.node->primitiveSet = ns.triangles;
			for (SampleInd& sample : ns.samples) {
				if (sample.prim >= 0) ns.node->insert(sample.prim, sampler->getSample(sample.raynr));
			}
			auto end_time1 = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1).count();
			std::cout << "Filled leaf in " << diff << " ms" << std::endl;
		}
	}

	static Split getRandomSplitterFromTriangles(RaySpaceTree* rst, std::set<int> tris) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r * (float)(tris.size()-1));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);
		std::set<int>::iterator it = tris.begin();
		std::advance(it, rtri);
		int triIndex = *it;
		int fromIndex = rst->model->indices[3 * triIndex + redge];
		int toIndex = rst->model->indices[3 * triIndex + (redge + 1) % 3];
		return { Ray(rst->model->vertices[fromIndex]->pos, rst->model->vertices[toIndex]->pos) };

	}

	static Split getRandomSplitterFromEdges(std::set<Edge*, Edge::cmp_ptr>& edges) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r * (float)(edges.size() - 1));
		std::set<Edge*>::iterator it = edges.begin();
		std::advance(it, rtri);
		Edge* e = *it;
		return { e->ray, e };
	}

	static std::vector<Split> getSplitterCandidates(RaySpaceTree* rst, int tries, NodeSamples& ns) {
		std::vector<Split> testRays = std::vector<Split>();

		if (ns.edges.size() <= tries) {
			for (Edge* e : ns.edges) 
				testRays.push_back({ e->ray, e });
			//for (int tri : tris) {
			//	for (int j = 0; j < 3; j++) {
			//		int fromIndex = rst->model->indices[3 * tri + j];
			//		int toIndex = rst->model->indices[3 * tri + (j + 1) % 3];
			//		Ray testRay = Ray(rst->model->verticesIndexed[fromIndex], rst->model->verticesIndexed[toIndex]);
			//		testRays.push_back(testRay);
			//	}
			//}
		}
		else {
			int i = 0; int j = 0;
			int upperlimit = 2 * tries;
			while (i < tries && j < upperlimit) {
				//Ray testRay = getRandomSplitterFromTriangles(rst, tris);
				Split testSplit = getRandomSplitterFromEdges(ns.edges);
				bool dupli = false;
				for (auto split : rst->splitters) if (split.ray.equal(testSplit.ray, 1E-6)) { dupli = true; break; }
				j++;
				if (dupli) continue;
				testRays.push_back(testSplit);
				i++;
			}
			if (testRays.size() == 0) std::cout << "Could not find any non-duplicate candidates" << std::endl;
		}
		return testRays;
	}


	static bool getBestSplitter(RaySpaceTree *rst, Split& splitter, Options::BuildOptions& options,
								NodeSamples& nsLeft, NodeSamples& nsRight, NodeSamples ns, 
								Sampler* sampler, int tries, bool print) {
		int t = (int)time(NULL);
		srand(t);

		float bestSAH = ns.triangles.size() * 1.f;
		if (print) std::cout << "Parent SAH = " << bestSAH << std::endl;
		float parentSAH = bestSAH;

		std::vector<Split> testSplitters = getSplitterCandidates(rst, tries, ns);

		// GET NODE AND THEN TRACE BACK WHICH SPLITTERS ARE ALREADY IN USE.
		// STORE IN SET OR IN VECTOR?
		// THEN ADD EVERYTIME SOMETHING NEW IS TRIED.
		// AND CHECK IF NEW OPTION IT'S ALREADY USED.

		for (Split& testSplitter : testSplitters) {
			NodeSamples nsLeftNew, nsRightNew;

			int pixL = 0;
			int pixR = 0;
			for (SampleInd sample : ns.samples) {
				Ray s_ray = sampler->getSample(sample.raynr);
				int tri = sample.prim;
				if (testSplitter.ray.side(s_ray)) {
					nsLeftNew.samples.push_back(sample);
					if (tri >= 0) {
						nsLeftNew.triangles.insert(tri);
						for (Edge* e : rst->model->triangles[tri]->edges) nsLeftNew.edges.insert(e);
					}
					pixL++;
				}
				else {
					nsRightNew.samples.push_back(sample);
					if (tri >= 0) {
						nsRightNew.triangles.insert(tri);
						for (Edge* e : rst->model->triangles[tri]->edges) nsRightNew.edges.insert(e);
					}
					pixR++;
				}
			}

			int samL = nsLeftNew.samples.size();
			int samR = nsRightNew.samples.size();

			float fracLeft = 1.f * (samL) / float(ns.samples.size());
			float fracRight = 1.f * (samR) / float(ns.samples.size());
			

			float sahleft = nsLeftNew.triangles.size() * fracLeft; // sah
			float sahright = nsRightNew.triangles.size() * fracRight;
			if (print) std::cout << "with " << fracLeft << ", " << fracRight << " samples (L,R) and " << 
									nsLeftNew.triangles.size() << ", " << nsRightNew.triangles.size() << 
									" triangles (L,R) for a SAH of " << sahleft + sahright;

			if (sahleft + sahright < bestSAH) {
				if (print) std::cout << " which is better!" << std::endl;
				bestSAH = sahleft + sahright;
				splitter = testSplitter;
				nsLeft = nsLeftNew;
				nsRight = nsRightNew;
			}
			else {
				if (print) std::cout << std::endl;
			}
		}
		if (print) std::cout << "Old value: " << ns.triangles.size() << ", New value: " << bestSAH << std::endl;
		if (ns.triangles.size() == bestSAH) return false;
		return true;
	}

};

