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
#include "shader.h"
#include "camera.h"
#include "rst.h"
#include "sample.h"

class RSTBuilderSamples : public RSTBuilder<RSTBuilderSamples> {
private:
	//static std::queue<nodeSamples> toProcess;
	//std::vector<Ray> raytraceRays;

public:

	static void build(Options::BuildOptions& options, RaySpaceTree* rst, VisComponents& visComp) {
		// spheresampler for now, will be cubesampling
		SphereSampler sampler = SphereSampler(&(rst->model->boundingSphere), options.noSamples);
		if (rst->alldir) sampler.createSamplesFull();
		else sampler.createSamples(rst->maindir);
		visComp.sampling_vao = sampler.vaoGeneration();
		visComp.sampling_size = sampler.samples.size();

		Orthocamera ocam(rst->model->radius, rst->model->radius, rst->model->radius * 2 > 30.f ? rst->model->radius * 2 : 30.f);
		Camera* cam = &ocam;
		TextureRenderer texrender;
		if (options.rasterizationSampling) {
			Shader rstProgram = Shader("rst.vert", "rst.frag");
			texrender = TextureRenderer(rstProgram, options.s_w, options.s_h);
		}

		if (options.construct == Options::ADAPTIVE) createAdaptive(rst, options, cam, sampler, texrender);
		else fill(rst, options, false, cam, sampler, texrender);

	}

	static glm::vec2 getPixPosFromIndex(int w, int h, int x, int y) {
		return { (x + .5f) / w * 2.0f - 1.0f, 1.0f - (y + .5f) / h * 2.0f };
	}

	static void makeSample(RaySpaceTree* rst, Options::BuildOptions& options, Camera* cam, bool store, GLfloat* pixels = 0,
							std::vector<SampleInd>& samples = std::vector<SampleInd>(), int camNum = 0,
							std::set<int>& triangles = std::set<int>()) {

		for (int y = 0; y < options.s_h; y++) {
			int yind = (options.s_h - 1 - y);
			for (int x = 0; x < options.s_w; x++) {
				Ray ray = cam->pixRayDirection(getPixPosFromIndex(options.s_w, options.s_h,x,y));
				if (!rst->alldir && !rst->model->boundingCube.intersectSide(rst->maindir, ray)) continue;

				int pixloc = yind * options.s_w + x;
				int tri = -1;
				if (options.rasterizationSampling) tri = int(pixels[pixloc]) - 1 ;
				else {
					float t = 0.f;
					rst->model->getIntersectionEmbree(ray, tri, t, true);
				}
				if (options.construct == Options::ADAPTIVE) {
					samples.push_back(SampleInd({ tri, camNum * options.s_h * options.s_w + pixloc }));
					if (tri >= 0) triangles.insert(tri);
				}
				else if (tri >= 0) rst->putPrimitive(ray, tri, store);
					// Extra checks, only necessary because rasterization not precise enough
					//if (rst->model->getIntersectionWithPrim(tri_id, ray, t))
					//else if (rst->model->getIntersectionEmbree(ray, tri_id, t))

			}
		}
	}

	static void fill(RaySpaceTree* rst, Options::BuildOptions& options, bool fillmore, Camera* cam, SphereSampler& sampler, TextureRenderer& texrender = TextureRenderer()) {

		std::cout << "Creating camera samples to fill tree..." << std::endl;
		int i = 0;
		while (i < sampler.samples.size()) {
			if (fillmore && i % sampler.ratio == 0) {
				i++;
				continue;
			}
			cam->setPositionAndForward(sampler.samples[i], rst->model->center);
			GLfloat* pixels = 0;
			if (options.rasterizationSampling) pixels = texrender.render(cam, rst->model);
			makeSample(rst, options, cam, !fillmore, pixels);
			i += sampler.ratio;
		}
	}

	static void createAdaptive(RaySpaceTree* rst, Options::BuildOptions& options, Camera* cam, SphereSampler& sampler, TextureRenderer& texrender = TextureRenderer()) {
		
		std::vector<SampleInd>& samples = std::vector<SampleInd>();
		std::set<int> triangles = std::set<int>();
		std::vector<Camera*> cams = std::vector<Camera*>();

		std::cout << "Rasterizing camera samples to create tree..." << std::endl;
		std::vector<Ray> raytraceRays;
		// prepare samples
		for (int i = 0; i < sampler.samples.size(); i += sampler.ratio) {
			Camera* cam1 = cam->makeCopy();
			cam1->setPositionAndForward(sampler.samples[i], rst->model->center);
			GLfloat* pixels = 0;
			if (options.rasterizationSampling) pixels = texrender.render(cam1, rst->model);
			makeSample(rst, options, cam1, options.storeRays, pixels, samples, i / sampler.ratio, triangles);
			cams.push_back(cam1);
		}
		std::cout << "Constructing the tree..." << std::endl;
		constructAdaptive(rst, options, cams, samples, triangles);
		//if (option == 0) 
		//else if (option == 1) rst->constructSmartRandom(texrender.res, cams, samples, tris);
	}

	static Ray getRayFromNum(std::vector<Camera*>& cams, int raynr,int w, int h) {
		int camNr = raynr / (w * h);
		int rayInCam = raynr % (w * h);
		int y = rayInCam / w;
		int x = rayInCam % w;
		const glm::vec2 pixpos = getPixPosFromIndex(w, h, x, y);
		return cams[camNr]->pixRayDirection(pixpos);
	}

	static void constructAdaptive(RaySpaceTree* rst, Options::BuildOptions& options, std::vector<Camera*>& cams, std::vector<SampleInd>& samples, std::set<int>& tris) {
		std::set<Edge*, Edge::cmp_ptr> edges;
		for (int i : tris) for (Edge* e : rst->model->triangles[i]->edges) edges.insert(e);


		std::queue<nodeSamples> toProcess;
		constructAdaptive(rst, 0, rst->rootNode, options, cams, samples, tris, edges, toProcess);
		samples = std::vector <SampleInd>();
		tris = std::set<int>();

		while (!toProcess.empty()) {
			nodeSamples ns = toProcess.front();
			constructAdaptive(rst, ns.level, ns.node, options, cams, ns.samples, ns.triangles, ns.edges, toProcess);
			toProcess.pop();

		}
	}

	static void constructAdaptive(RaySpaceTree *rst, int level, Node* node, Options::BuildOptions& options, std::vector<Camera*>& cams,
									std::vector<SampleInd>& samples, std::set<int>& tris, std::set<Edge*, Edge::cmp_ptr>& edges, std::queue<nodeSamples>& toProcess) {
		auto start_time1 = std::chrono::high_resolution_clock::now();

		std::cout << "level: " << level << std::endl;
		bool checkmore = false;
		if (tris.size() == 1) {
			for (auto s : samples)
				if (s.prim < 0) {
					checkmore = true;
					break;
				}
		}

		if (level <= rst->depth, (tris.size() > 1 || checkmore)) {
			std::cout << "Node with " << samples.size() << " samples" << std::endl;
			std::cout << "Finding best splitter..." << std::endl;

			auto start_time = std::chrono::high_resolution_clock::now();
			Ray splitter;

			//change this to all nodesamples
			std::vector<SampleInd> samplesL, samplesR;
			std::set<int> trisL, trisR; 
			std::set<Edge*, Edge::cmp_ptr> edgesL, edgesR;

			if (getBestSplitter(rst, splitter, options, samplesL, trisL, edgesL, samplesR, trisR, edgesR, samples, tris, edges, cams, 10,true)) {

				auto end_time = std::chrono::high_resolution_clock::now();
				auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
				std::cout << "Found splitter in " << diff << " ms" << std::endl;

				node->splitter = splitter;
				rst->splitters.push_back(splitter);

				Node* leftnode = new Node(node->index * 2 + 1, level + 1);
				node->leftNode = leftnode;
				leftnode->parent = node;
				rst->nodes.push_back(leftnode);
				nodeSamples nsL = { trisL, samplesL, edgesL, node->leftNode, level + 1 };
				toProcess.push(nsL);
				//constructAdaptive(level + 1, node->leftNode, res, cams, samplesL, trisL, print);

				Node* rightnode = new Node(node->index * 2 + 2, level + 1);
				rightnode->parent = node;
				node->rightNode = rightnode;
				rst->nodes.push_back(rightnode);
				nodeSamples nsR = { trisR, samplesR, edgesL, node->rightNode, level + 1 };
				toProcess.push(nsR);
				//constructAdaptive(level + 1, node->rightNode, res, cams, samplesR, trisR, print);
				return;
			}
		}
		std::cout << "Leaf with " << samples.size() << " samples and " << tris.size() << " triangles" << std::endl;
		
		rst->noLeaves++;
		node->leaf = true;

		if (tris.size() > 0) {
			node->primitiveSet = tris;
			for (SampleInd& sample : samples) {
				if (sample.prim >= 0) node->insert(sample.prim, getRayFromNum(cams, sample.raynr, options.s_w, options.s_h));
			}
			auto end_time1 = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1).count();
			std::cout << "Filled leaf in " << diff << " ms" << std::endl;
		}
	}

	static Ray getRandomSplitterFromTriangles(RaySpaceTree* rst, std::set<int> tris) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r * (float)(tris.size()-1));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);
		std::set<int>::iterator it = tris.begin();
		std::advance(it, rtri);
		int triIndex = *it;
		int fromIndex = rst->model->indices[3 * triIndex + redge];
		int toIndex = rst->model->indices[3 * triIndex + (redge + 1) % 3];
		return Ray(rst->model->vertices[fromIndex]->pos, rst->model->vertices[toIndex]->pos);

	}

	static Ray getRandomSplitterFromEdges(std::set<Edge*, Edge::cmp_ptr>& edges) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r * (float)(edges.size() - 1));
		std::set<Edge*>::iterator it = edges.begin();
		std::advance(it, rtri);
		Edge* e = *it;
		return e->ray;
	}

	static std::vector<Ray> getSplitterCandidates(RaySpaceTree* rst, int tries, std::set<int>& tris, std::set<Edge*, Edge::cmp_ptr>& edges) {
		std::vector<Ray> testRays = std::vector<Ray>();

		if (edges.size() <= tries) {
			for (Edge* e : edges) testRays.push_back(e->ray);
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
				Ray testRay = getRandomSplitterFromEdges(edges);
				bool dupli = false;
				for (auto r : rst->splitters) if (r.equal(testRay, 1E-6)) { dupli = true; break; }
				j++;
				if (dupli) continue;
				testRays.push_back(testRay);
				i++;
			}
			if (testRays.size() == 0) std::cout << "Could not find any non-duplicate candidates" << std::endl;
		}
		return testRays;
	}


	static bool getBestSplitter(RaySpaceTree *rst, Ray& splitter, Options::BuildOptions& options,
								std::vector<SampleInd>& samplesLx, std::set<int>& trisLx, std::set<Edge*, Edge::cmp_ptr>& edgesLx,
								std::vector<SampleInd>& samplesRx, std::set<int>& trisRx, std::set<Edge*, Edge::cmp_ptr>& edgesRx,
								std::vector<SampleInd>& samples, std::set<int>& tris, std::set<Edge*, Edge::cmp_ptr>& edges,
								std::vector<Camera*>& cams, int tries, bool print) {
		int t = (int)time(NULL);
		srand(t);

		float bestSAH = tris.size() * 1.f;
		if (print) std::cout << "Parent SAH = " << bestSAH << std::endl;
		float parentSAH = bestSAH;

		std::vector<Ray> testSplitters = getSplitterCandidates(rst, tries, tris, edges);

		// GET NODE AND THEN TRACE BACK WHICH SPLITTERS ARE ALREADY IN USE.
		// STORE IN SET OR IN VECTOR?
		// THEN ADD EVERYTIME SOMETHING NEW IS TRIED.
		// AND CHECK IF NEW OPTION IT'S ALREADY USED.

		for (Ray& testSplitter : testSplitters) {
			std::vector<SampleInd> samplesL, samplesR;
			std::set<int> trisL, trisR;
			std::set<Edge*, Edge::cmp_ptr> edgesL, edgesR;

			int pixL = 0;
			int pixR = 0;
			for (SampleInd sample : samples) {
				Ray s_ray = getRayFromNum(cams, sample.raynr, options.s_w, options.s_h);
				int tri = sample.prim;
				if (testSplitter.side(s_ray)) {
					samplesL.push_back(sample);
					if (tri >= 0) {
						trisL.insert(tri);
						for (Edge* e : rst->model->triangles[tri]->edges) edgesL.insert(e);
					}
					pixL++;
				}
				else {
					samplesR.push_back(sample);
					if (tri >= 0) {
						trisR.insert(tri);
						for (Edge* e : rst->model->triangles[tri]->edges) edgesL.insert(e);
					}
					pixR++;
				}
			}

			int samL = samplesL.size();
			int samR = samplesR.size();

			float fracLeft = 1.f * (samL) / float(samples.size());
			float fracRight = 1.f * (samR) / float(samples.size());
			

			float sahleft = trisL.size() * fracLeft; // sah
			float sahright = trisR.size() * fracRight;
			if (print) std::cout << "with " << fracLeft << ", " << fracRight << " samples (L,R) and " << trisL.size() << ", " << trisR.size() << " triangles (L,R) for a SAH of " << sahleft + sahright;

			if (sahleft + sahright < bestSAH) {
				if (print) std::cout << " which is better!" << std::endl;
				bestSAH = sahleft + sahright;
				splitter = testSplitter;
				samplesLx = samplesL;
				samplesRx = samplesR;
				trisLx = trisL;
				trisRx = trisR;
				edgesLx = edgesL;
				edgesRx = edgesR;
			}
			else {
				if (print) std::cout << std::endl;
			}
		}
		if (print) std::cout << "Old value: " << tris.size() << ", New value: " << bestSAH << std::endl;
		if (tris.size() == bestSAH) return false;
		return true;
	}


	//void constructAdaptive(glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples,
	//	std::set<int>& tris, bool print);
	//void constructAdaptive(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples,
	//	std::set<int>& tris, bool print, std::vector<Ray> splitters);
	//void constructAdaptive2(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples,
	//	std::set<int>& tris, int totSampSize);
	//
	//void constructSmartRandom(glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples, std::set<int>& tris);
	//void constructSmartRandom(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples,
	//	std::set<int>& tris, std::vector<Ray> splitters);
	//
	//bool getBestSplitter(Ray& splitter, glm::ivec2 res, std::vector<std::pair<int, int>>& samplesL, std::set<int>& trisL,
	//	std::vector<std::pair<int, int>>& samplesR, std::set<int>& trisR, std::vector<Camera*>& cams, int tries,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray>& splitters);
	//Ray getBestSplitter2(glm::ivec2 res, std::vector<Camera*>& cams, int tries,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize);
	//
	//void RaySpaceTree::constructSmartRandom(glm::ivec2 res, std::vector<Camera*>& cams,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris) {
	//	int t = (int)time(NULL);
	//	srand(t);
	//	std::vector<Ray> splitters;
	//	constructSmartRandom(0, rootNode, res, cams, samples, tris, splitters);
	//}

	//void RaySpaceTree::constructSmartRandom(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, std::vector<Ray> splitters) {
	//	if (level >= depth || tris.size() <= 1) {
	//		node->leaf = true;
	//		noLeaves++;
	//		if (tris.size() > 0) {
	//			node->primitiveSet = tris;
	//			for (auto sample : samples) {
	//				if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
	//			}
	//		}
	//		return;
	//	}
	//	else {
	//		float r = (float)rand() / static_cast <float> (RAND_MAX);
	//		int rtri = int(r * (float)tris.size());
	//		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
	//		int redge = int(r2 * 3.f);
	//		std::set<int>::iterator it = tris.begin();
	//		std::advance(it, rtri);
	//		int triIndex = *it;
	//		int fromIndex = model->indices[3 * triIndex + redge];
	//		int toIndex = model->indices[3 * triIndex + (redge + 1) % 3];
	//		Ray splitter = Ray(model->verticesIndexed[fromIndex], model->verticesIndexed[toIndex]);
	//
	//		bool dupli = false;
	//		for (Ray r : splitters) if (splitter.equal(r, 1E-6)) { dupli = true; break; }
	//		if (dupli) constructSmartRandom(level, node, res, cams, samples, tris, splitters);
	//
	//		node->splitter = splitter;
	//		splitters.push_back(splitter);
	//
	//		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
	//		std::set<int> trisL = std::set<int>();
	//		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
	//		std::set<int> trisR = std::set<int>();
	//
	//		for (auto sample : samples) {
	//			bool left = intoLeftNode(splitter, cams, sample.first, res);
	//			int tri = sample.second;
	//			if (left) {
	//				samplesL.push_back(sample);
	//				if (tri >= 0) trisL.insert(tri);
	//			}
	//			else {
	//				samplesR.push_back(sample);
	//				if (tri >= 0) trisR.insert(tri);
	//			}
	//		}
	//
	//		Node* leftnode = new Node(node->index * 2 + 1, level + 1);
	//		node->leftNode = leftnode;
	//		leftnode->parent = node;
	//		nodes.push_back(leftnode);
	//		constructSmartRandom(level + 1, node->leftNode, res, cams, samplesL, trisL, splitters);
//
	//		Node* rightnode = new Node(node->index * 2 + 2, level + 1);
	//		rightnode->parent = node;
	//		node->rightNode = rightnode;
	//		nodes.push_back(rightnode);
	//		constructSmartRandom(level + 1, node->rightNode, res, cams, samplesR, trisR, splitters);
	//	}
	//}

	//void RaySpaceTree::constructAdaptive(glm::ivec2 res, std::vector<Camera*>& cams,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print)
	//{
	//	std::vector<Ray> splitters;
	//
	//	constructAdaptive(0, rootNode, res, cams, samples, tris, print, splitters);
	//	samples = std::vector < std::pair<int, int>>();
	//	tris = std::set<int>();
	//
	//	while (!toProcess.empty()) {
	//		nodeSamples ns = toProcess.front();
	//		constructAdaptive(ns.level, ns.node, res, cams, ns.samples, ns.triangles, print, splitters);
	//		toProcess.pop();
	//	}
	//}

	//void RaySpaceTree::constructAdaptive2(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize) {
	//	auto start_time = std::chrono::high_resolution_clock::now();
	//	std::cout << "level: " << level << std::endl;
	//	if (tris.size() > 1) {
	//		std::cout << "Node with " << samples.size() << " samples and " << tris.size() << " triangles." << std::endl;
	//		std::cout << "Finding best splitter..." << std::endl;
	//
	//		auto start_time = std::chrono::high_resolution_clock::now();
	//		Ray splitter = getBestSplitter2(res, cams, 10, samples, tris, totSampSize);
	//		auto end_time = std::chrono::high_resolution_clock::now();
	//		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	//		std::cout << "Found splitter in " << diff << " ms" << std::endl;
	//
	//		node->splitter = splitter;
	//		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
	//		std::set<int> trisL = std::set<int>();
	//		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
	//		std::set<int> trisR = std::set<int>();
	//
	//		std::cout << "Putting samples left and right..." << std::endl;
	//		start_time = std::chrono::high_resolution_clock::now();
	//		for (auto sample : samples) {
	//			bool left = intoLeftNode(splitter, cams, sample.first, res);
	//			int tri = sample.second;
	//			if (left) {
	//				samplesL.push_back(sample);
	//				if (tri >= 0) trisL.insert(tri);
	//			}
	//			else {
	//				samplesR.push_back(sample);
	//				if (tri >= 0) trisR.insert(tri);
	//			}
	//		}
	//		end_time = std::chrono::high_resolution_clock::now();
	//		diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	//		std::cout << "Samples put left and right in " << diff << " ms" << std::endl;
	//
	//		Node* leftnode = new Node(node->index * 2 + 1, level + 1);
	//		node->leftNode = leftnode;
	//		leftnode->parent = node;
	//		nodes.push_back(leftnode);
	//		constructAdaptive2(level + 1, node->leftNode, res, cams, samplesL, trisL, totSampSize);
	//
	//		Node* rightnode = new Node(node->index * 2 + 2, level + 1);
	//		rightnode->parent = node;
	//		node->rightNode = rightnode;
	//		nodes.push_back(rightnode);
	//		constructAdaptive2(level + 1, node->rightNode, res, cams, samplesR, trisR, totSampSize);
	//
	//		return;
	//	}
	//	noLeaves++;
	//	std::cout << "Leaf with " << samples.size() << " samples and " << tris.size() << " triangles" << std::endl;
	//	node->leaf = true;
	//	if (tris.size() > 0) {
	//		node->primitiveSet = tris;
	//		for (auto sample : samples) {
	//			if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
	//		}
	//		auto end_time = std::chrono::high_resolution_clock::now();
	//		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	//		std::cout << "Filled leaf in " << diff << " ms" << std::endl;
	//	}
	//}

	//void RaySpaceTree::constructAdaptive(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray> splitters) {
	//	auto start_time = std::chrono::high_resolution_clock::now();
	//	if (print) std::cout << "level: " << level << std::endl;
	//	bool checkmore = false;
	//	if (tris.size() == 1) {
	//		for (auto s : samples)
	//			if (s.second < 0) checkmore = true;
	//	}
	//
	//	if (tris.size() > 1 || checkmore) {
	//		if (print) std::cout << "Node with " << samples.size() << " samples" << std::endl;
	//		if (print) std::cout << "Finding best splitter..." << std::endl;
	//		auto start_time = std::chrono::high_resolution_clock::now();
	//		Ray splitter;
	//		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
	//		std::set<int> trisL = std::set<int>();
	//		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
	//		std::set<int> trisR = std::set<int>();
	//		if (getBestSplitter(splitter, res, samplesL, trisL, samplesR, trisR, cams, 10, samples, tris, print, splitters)) {
	//			auto end_time = std::chrono::high_resolution_clock::now();
	//			auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	//			if (print) std::cout << "Found splitter in " << diff << " ms" << std::endl;
	//			node->splitter = splitter;
	//			splitters.push_back(splitter);
	//			Node* leftnode = new Node(node->index * 2 + 1, level + 1);
	//			node->leftNode = leftnode;
	//			leftnode->parent = node;
	//			nodes.push_back(leftnode);
	//			nodeSamples nsL = { trisL, samplesL, node->leftNode, level + 1 };
	//			toProcess.push(nsL);
	//			//constructAdaptive(level + 1, node->leftNode, res, cams, samplesL, trisL, print);
	//			Node* rightnode = new Node(node->index * 2 + 2, level + 1);
	//			rightnode->parent = node;
	//			node->rightNode = rightnode;
	//			nodes.push_back(rightnode);
	//			nodeSamples nsR = { trisR, samplesR, node->rightNode, level + 1 };
	//			toProcess.push(nsR);
	//			//constructAdaptive(level + 1, node->rightNode, res, cams, samplesR, trisR, print);
	//			return;
	//		}
	//	}
	//	noLeaves++;
	//	if (print) std::cout << "Leaf with " << samples.size() << " samples and " << tris.size() << " triangles" << std::endl;
	//	node->leaf = true;
	//	if (tris.size() > 0) {
	//		node->primitiveSet = tris;
	//		for (auto sample : samples) {
	//			if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
	//		}
	//		auto end_time = std::chrono::high_resolution_clock::now();
	//		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	//		if (print) std::cout << "Filled leaf in " << diff << " ms" << std::endl;
	//	}
	//}

	//Ray RaySpaceTree::getBestSplitter2(glm::ivec2 res, std::vector<Camera*>& cams, int tries,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize) {
	//	int t = (int)time(NULL);
	//	srand(t);
	//	float bestTriFrac = 0;
	//	Ray splitter;
	//
	//	//for (int i = 0; i < tris.size(); i++) {
	//		//for (int j = 0; j < 3; j++) {
	//	for (int i = 0; i < tries; i++) {
	//		std::cout << "Checking line no " << i << "... ";
	//
	//		// choose the first edge of a random triangle in this node
	//		float r = (float)rand() / static_cast <float> (RAND_MAX);
	//		int rtri = int(r * (float)tris.size());
	//		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
	//		int redge = int(r2 * 3.f);
	//		std::set<int>::iterator it = tris.begin();
	//		std::advance(it, rtri);
	//		int tri = *it;
	//		Ray testRay = Ray(model->vertices[3 * tri + redge].pos, model->vertices[3 * tri + (redge + 1) % 3].pos);
	//		//Ray testRay = Ray(vertices[3 * i + j].pos, vertices[3 * i + (j + 1) % 3].pos);
	//
	//
	//		std::set<int> trisLeft = std::set<int>();
	//		std::set<int> trisRight = std::set<int>();
	//		int zeroL = 0; int zeroR = 0;
	//
	//		for (auto sample : samples) {
	//			bool left = intoLeftNode(testRay, cams, sample.first, res);
	//			int tri = sample.second;
	//			if (left) {
	//				if (tri >= 0) {
	//					trisLeft.insert(tri);
	//
	//				}
	//				else zeroL++;
	//			}
	//			else {
	//				if (tri >= 0) {
	//					trisRight.insert(tri);
	//				}
	//				else zeroR++;
	//			}
	//		}
	//
	//		float frac = 0;
	//		if (trisLeft.size() && trisRight.size()) frac = 1.f * trisLeft.size() / (1.f * trisRight.size());
	//		if (frac > 1.f) frac = 1.f / frac;
	//		std::cout << " for a tri frac of " << frac;
	//
	//		if (frac > bestTriFrac) {
	//			std::cout << " which is better!" << std::endl;
	//			bestTriFrac = frac;
	//			splitter = testRay;
	//		}
	//		else {
	//			std::cout << std::endl;
	//		}
	//		//	}
	//	}
	//	std::cout << "New value: " << bestTriFrac << std::endl;
	//	//if (bestTriFrac == 0) getBestSplitter(splitter, res, cams, tries, samples, tris, totSampSize);
	//	return splitter;
	//}

	//bool RaySpaceTree::getBestSplitter(Ray& splitter, glm::ivec2 res, std::vector<std::pair<int, int>>& samplesLx, std::set<int>& trisLx,
	//	std::vector<std::pair<int, int>>& samplesRx, std::set<int>& trisRx, std::vector<Camera*>& cams, int tries,
	//	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray>& splitters) {
	//	int t = (int)time(NULL);
	//	srand(t);
	//
	//	float bestSAH = tris.size() * 1.f;
	//	if (print) std::cout << "Parent SAH = " << bestSAH << std::endl;
	//	float parentSAH = bestSAH;
	//
	//	std::vector<Ray> testRays = std::vector<Ray>();
	//	std::vector<int> triIndices = std::vector<int>();
	//	std::vector<int> triIndices2 = std::vector<int>();
	//
	//	if (tris.size() <= 3) {
	//		for (int tri : tris) {
	//			for (int j = 0; j < 3; j++) {
	//				int triIndex = tri;
	//				int fromIndex = model->indices[3 * tri + j];
	//				int toIndex = model->indices[3 * tri + (j + 1) % 3];
	//				Ray testRay = Ray(model->verticesIndexed[fromIndex], model->verticesIndexed[toIndex]);
	//
	//				//int triIndex2 = -5;
	//				//if (model->triPerVertex[fromIndex].size() > 1 && model->triPerVertex[toIndex].size() > 1) {
	//				//	for (auto trix : model->triPerVertex[fromIndex]) {
	//				//		if (triIndex2 >= 0) break;
	//				//		for (auto triy : model->triPerVertex[fromIndex]) {
	//				//			if (trix != triIndex && trix == triy) {
	//				//				triIndex2 = trix;
	//				//				break;
	//				//			}
	//				//		}
	//				//	}
	//				//}
	//				testRays.push_back(testRay);
	//				//triIndices.push_back(triIndex);
	//				//.push_back(triIndex2);
	//			}
	//		}
	//	}
	//	else {
	//		int i = 0;
	//		int upperlimit = 2 * tries;
	//		while (i < tries && i < upperlimit) {
	//			float r = (float)rand() / static_cast <float> (RAND_MAX);
	//			int rtri = int(r * (float)tris.size());
	//			float r2 = (float)rand() / static_cast <float> (RAND_MAX);
	//			int redge = int(r2 * 3.f);
	//			std::set<int>::iterator it = tris.begin();
	//			std::advance(it, rtri);
	//			int triIndex = *it;
	//			int fromIndex = model->indices[3 * triIndex + redge];
	//			int toIndex = model->indices[3 * triIndex + (redge + 1) % 3];
	//			Ray testRay = Ray(model->verticesIndexed[fromIndex], model->verticesIndexed[toIndex]);
	//			bool dupli = false;
	//			for (auto r : splitters) if (r.equal(testRay, 1E-6)) { dupli = true; break; }
	//			if (dupli) continue;
	//			i++;
	//			testRays.push_back(testRay);
	//
	//			//int triIndex2 = -5;
	//			//if (model->triPerVertex[fromIndex].size() > 1 && model->triPerVertex[toIndex].size() > 1) {
	//			//	for (auto trix : model->triPerVertex[fromIndex]) {
	//			//		if (triIndex2 >= 0) break;
	//			//		for (auto triy : model->triPerVertex[fromIndex]) {
	//			//			if (trix != triIndex && trix == triy) {
	//			//				triIndex2 = trix;
	//			//				break;
	//			//			}
	//			//		}
	//			//	}
	//			//}	
	//			//triIndices.push_back(triIndex);
	//			//triIndices2.push_back(triIndex2);
	//		}
	//		if (testRays.size() == 0) std::cout << "Could not find any non-duplicate candidates" << std::endl;
	//	}
	//
	//	// GET NODE AND THEN TRACE BACK WHICH SPLITTERS ARE ALREADY IN USE.
	//	// STORE IN SET OR IN VECTOR?
	//	// THEN ADD EVERYTIME SOMETHING NEW IS TRIED.
	//	// AND CHECK IF NEW OPTION IT'S ALREADY USED.
	//
	//	// What if one of the samples incorrectly hits empty space? -- with the additional raytracing that should
	//	// now be the only check that needs to be done.
	//	for (int i = 0; i < testRays.size(); i++) {
	//		//std::vector<int> checktri = { triIndices[i], triIndices2[i] };
	//		//std::vector<std::vector<int>> checkIndTris = std::vector<std::vector<int>>(2*checktri.size());
	//
	//		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
	//		std::set<int> trisL = std::set<int>();
	//		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
	//		std::set<int> trisR = std::set<int>();
	//
	//		int pixL = 0;
	//		int pixR = 0;
	//		for (auto sample : samples) {
	//			bool left = intoLeftNode(testRays[i], cams, sample.first, res);
	//			int tri = sample.second;
	//
	//			//for (int q = 0; q < checktri.size(); q++) {
	//			//	if (checktri[q] == tri) {
	//			//		if (left) checkIndTris[2*q].push_back(pixL);
	//			//		else checkIndTris[2*q+1].push_back(pixR);;
	//			//	}
	//			//}
	//			if (left) {
	//				samplesL.push_back(sample);
	//				if (tri >= 0) trisL.insert(tri);
	//				pixL++;
	//			}
	//			else {
	//				samplesR.push_back(sample);
	//				if (tri >= 0) trisR.insert(tri);
	//				pixR++;
	//			}
	//		}
	//
	//		int samL = samplesL.size();
	//		int samR = samplesR.size();
	//		//std::vector<int> toRemoveL;
	//		//std::vector<int> toRemoveR;
	//
	//		//for (int q = 0; q < checktri.size(); q++) {
	//		//	int leftsize = checkIndTris[2 * q].size();
	//		//	int rightsize = checkIndTris[2 * q + 1].size();
	//		//	if (leftsize > 0 && rightsize > 0) {
	//		//		if (leftsize < rightsize) {
	//		//			//for (auto sp : checkIndTris[2 * q]) {
	//		//			//	samplesR.push_back(samplesL[sp]);
	//		//			//	samplesL[sp] = samplesL.back();
	//		//			//	//samplesL.pop_back();
	//		//			//}
	//		//			samL--;
	//		//			samR++;
	//		//			trisL.erase(checktri[q]);
	//		//		}
	//		//		else {
	//		//			//for (auto sp : checkIndTris[2 * q + 1]) {
	//		//				//samplesL.push_back(samplesR[sp]);
	//		//				//samplesR[sp] = samplesR.back();
	//		//				//samplesR.pop_back();
	//		//			//}
	//		//			samL++;
	//		//			samR--;
	//		//			trisR.erase(checktri[q]);
	//		//		}
	//		//	}
	//		//}
	//
	//		float fracLeft = 1.f * (samL) / float(samples.size());
	//		float fracRight = 1.f * (samR) / float(samples.size());
	//
	//		float sahleft = trisL.size() * fracLeft; // sah
	//		float sahright = trisR.size() * fracRight;
	//		if (print) std::cout << "with " << fracLeft << ", " << fracRight << " samples (L,R) and " << trisL.size() << ", " << trisR.size() << " triangles (L,R) for a SAH of " << sahleft + sahright;
	//
	//		if (sahleft + sahright < bestSAH) {
	//			if (print) std::cout << " which is better!" << std::endl;
	//			bestSAH = sahleft + sahright;
	//			splitter = testRays[i];
	//			samplesLx = samplesL;
	//			samplesRx = samplesR;
	//			trisLx = trisL;
	//			trisRx = trisR;
	//		}
	//		else {
	//			if (print) std::cout << std::endl;
	//		}
	//	}
	//	if (print) std::cout << "Old value: " << tris.size() << ", New value: " << bestSAH << std::endl;
	//	if (tris.size() == bestSAH) return false;
	//	return true;
	//}
};

