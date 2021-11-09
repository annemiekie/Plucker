#include "rst.h"
#include <chrono>

void RaySpaceTree::construct(int option, std::vector<Ray>& rays) {
	int t = (int)time(NULL);
	//srand(t);
	std::vector<Ray> splitters;
	construct(0, rootNode, option, rays, splitters);
	noLeaves = pow((int)2, depth);
}

void RaySpaceTree::construct(int lvl, Node* node, int option, std::vector<Ray>& rays, std::vector<Ray> splitters) {
	if (lvl >= depth) {
		node->leaf = true;
		return;
	}
	Ray splitter;
	if (option == 0) { // edges
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * (model->vertices.size() / 3));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);
		splitter = Ray(model->vertices[3 * r1 + redge].pos, model->vertices[3 * r1 + (redge + 1) % 3].pos);
	}
	else if (option == 1) { // orthogonal from vertex
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * model->vertices.size());
		glm::vec3 pos = model->vertices[r1].pos;
		if (rand() % 2) splitter = Ray(glm::vec3(pos.x, 0, pos.z), glm::vec3(pos.x, 1, pos.z));
		else splitter = Ray(glm::vec3(pos.x, pos.y, 0), glm::vec3(pos.x, pos.y, 1));
	}
	else if (option == 2) { // orthogonal random
		float rX = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		float rZ = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		float rY = 2.f * float(rand()) / static_cast <float> (RAND_MAX);

		if (rand() % 2) splitter = Ray(glm::vec3(rX, 0, rZ), glm::vec3(rX, 1, rZ));
		else splitter = Ray(glm::vec3(rX, rY, 0), glm::vec3(rX, rY, 1));
	}
	else if (option == 3) {
		splitter = rays[lvl];
	}
	else if (option == 4) { // random offset to edges
		float r1 = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r1 * (model->vertices.size() / 3));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);
		float r3 = (float)rand() / static_cast <float> (RAND_MAX);
		float roffset = 0.1 * r3;
		splitter = Ray(model->vertices[3 * rtri + redge].pos + roffset, model->vertices[3 * rtri + (redge + 1) % 3].pos);
	}
	//else if (option == 5) { // random vertices?
		//float r = (float)rand() / static_cast <float> (RAND_MAX);
		//int r1 = int(r * (model->vertices.size()));
		//float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		//int r2 = int(r * (model->vertices.size()));
		//float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		//int redge = int(r2 * 3.f);
		//float r3 = (float)rand() / static_cast <float> (RAND_MAX);
		//float roffset = 0.1 * r3;
		//splitter = Ray(model->vertices[3 * r1 + redge].pos + roffset, model->vertices[3 * r1 + (redge + 1) % 3].pos);
	//}
	bool dupli = false;
	for (Ray r : splitters) if (splitter.equal(r, 1E-6)) {dupli = true; break;}
	if (dupli) construct(lvl, node, option, rays, splitters);
	else {
		splitters.push_back(splitter);

		node->splitter = splitter;
		Node* leftnode = new Node(node->index * 2 + 1, lvl + 1);
		node->leftNode = leftnode;
		leftnode->parent = node;
		nodes.push_back(leftnode);
		construct(lvl + 1, node->leftNode, option, rays, splitters);

		Node* rightnode = new Node(node->index * 2 + 2, lvl + 1);
		rightnode->parent = node;
		node->rightNode = rightnode;
		nodes.push_back(rightnode);
		construct(lvl + 1, node->rightNode, option, rays, splitters);
	}
}

void RaySpaceTree::fillExact() {
	std::vector<Ray> splitters;
	std::set<int> triangles;
	Ray r;
	if (alldir) {
		for (Node* node : nodes) {
			if (node->leaf) {
				for (int i = 0; i < model->primsize; i++) {
					if (check1Prim(i, r, node, false, 0)) node->primitiveSet.insert(i);
				}
			}
		}
	}
	else {
		#pragma omp parallel
		for (int i = 0; i < model->primsize; i++) {
			//if (!alldir) 
			if (check1Prim(i, r, rootNode, false, 0)) rootNode->primitiveSet.insert(i);//, 0
			//else  rootNode->primitiveSet.insert(i);
		}
		fillExact(rootNode->leftNode);
		fillExact(rootNode->rightNode);
	}

}

void RaySpaceTree::fillExact(Node* node) {
	Ray r;
	#pragma omp parallel
	for (int i: node->parent->primitiveSet) if (check1Prim(i, r, node, false, 1)) node->primitiveSet.insert(i);//, 1

	if (!node->leaf) {
		fillExact(node->leftNode);
		fillExact(node->rightNode);
	}
}

void RaySpaceTree::constructSmartRandom(glm::ivec2 res, std::vector<Camera*>& cams,
										std::vector<std::pair<int, int>>& samples, std::set<int>& tris) {
	int t = (int)time(NULL);
	srand(t);
	std::vector<Ray> splitters;
	constructSmartRandom(0, rootNode, res, cams, samples, tris, splitters);
}

void RaySpaceTree::constructSmartRandom(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, std::vector<Ray> splitters) {
	if (level >= depth || tris.size() <= 1){
		node->leaf = true;
		noLeaves++;
		if (tris.size() > 0) {
			node->primitiveSet = tris;
			for (auto sample : samples) {
				if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
			}
		}
		return;
	}
	else {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r * (float)tris.size());
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);
		std::set<int>::iterator it = tris.begin();
		std::advance(it, rtri);
		int triIndex = *it;
		int fromIndex = model->indices[3 * triIndex + redge];
		int toIndex = model->indices[3 * triIndex + (redge + 1) % 3];
		Ray splitter = Ray(model->vertices2[fromIndex], model->vertices2[toIndex]);

		bool dupli = false;
		for (Ray r : splitters) if (splitter.equal(r, 1E-6)) { dupli = true; break; }
		if (dupli) constructSmartRandom(level, node, res, cams, samples, tris, splitters);

		node->splitter = splitter;
		splitters.push_back(splitter);

		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
		std::set<int> trisL = std::set<int>();
		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
		std::set<int> trisR = std::set<int>();

		for (auto sample : samples) {
			bool left = intoLeftNode(splitter, cams, sample.first, res);
			int tri = sample.second;
			if (left) {
				samplesL.push_back(sample);
				if (tri >= 0) trisL.insert(tri);
			}
			else {
				samplesR.push_back(sample);
				if (tri >= 0) trisR.insert(tri);
			}
		}

		Node* leftnode = new Node(node->index * 2 + 1, level + 1);
		node->leftNode = leftnode;
		leftnode->parent = node;
		nodes.push_back(leftnode);
		constructSmartRandom(level + 1, node->leftNode, res, cams, samplesL, trisL, splitters);

		Node* rightnode = new Node(node->index * 2 + 2, level + 1);
		rightnode->parent = node;
		node->rightNode = rightnode;
		nodes.push_back(rightnode);
		constructSmartRandom(level + 1, node->rightNode, res, cams, samplesR, trisR, splitters);
	}
}

void RaySpaceTree::constructAdaptive(glm::ivec2 res, std::vector<Camera*>& cams,
	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print)
{
	std::vector<Ray> splitters;

	constructAdaptive(0, rootNode, res, cams, samples, tris, print, splitters);
	samples = std::vector < std::pair<int, int>>();
	tris = std::set<int>();

	while (!toProcess.empty()) {
		nodeSamples ns = toProcess.front();
		constructAdaptive(ns.level, ns.node, res, cams, ns.samples, ns.triangles, print, splitters);
		toProcess.pop();
	}
}

void RaySpaceTree::constructAdaptive2(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
									  std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize) {
	auto start_time = std::chrono::high_resolution_clock::now();
	std::cout << "level: " << level << std::endl;
	if (tris.size() > 1) {
		std::cout << "Node with " << samples.size() << " samples and " << tris.size() << " triangles." <<std::endl;
		std::cout << "Finding best splitter..." << std::endl;

		auto start_time = std::chrono::high_resolution_clock::now();
		Ray splitter = getBestSplitter2(res, cams, 10, samples, tris, totSampSize);
		auto end_time = std::chrono::high_resolution_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
		std::cout << "Found splitter in " << diff << " ms" << std::endl;

		node->splitter = splitter;
		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
		std::set<int> trisL = std::set<int>();
		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
		std::set<int> trisR = std::set<int>();

		std::cout << "Putting samples left and right..." << std::endl;
		start_time = std::chrono::high_resolution_clock::now();
		for (auto sample : samples) {
			bool left = intoLeftNode(splitter, cams, sample.first, res);
			int tri = sample.second;
			if (left) {
				samplesL.push_back(sample);
				if (tri >= 0) trisL.insert(tri);
			}
			else {
				samplesR.push_back(sample);
				if (tri >= 0) trisR.insert(tri);
			}
		}
		end_time = std::chrono::high_resolution_clock::now();
		diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
		std::cout << "Samples put left and right in " << diff << " ms" << std::endl;

		Node* leftnode = new Node(node->index * 2 + 1, level + 1);
		node->leftNode = leftnode;
		leftnode->parent = node;
		nodes.push_back(leftnode);
		constructAdaptive2(level + 1, node->leftNode, res, cams, samplesL, trisL, totSampSize);

		Node* rightnode = new Node(node->index * 2 + 2, level + 1);
		rightnode->parent = node;
		node->rightNode = rightnode;
		nodes.push_back(rightnode);
		constructAdaptive2(level + 1, node->rightNode, res, cams, samplesR, trisR, totSampSize);

		return;
	}
	noLeaves++;
	std::cout << "Leaf with " << samples.size() << " samples and " << tris.size() << " triangles" << std::endl;
	node->leaf = true;
	if (tris.size() > 0) {
		node->primitiveSet = tris;
		for (auto sample : samples) {
			if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
		std::cout << "Filled leaf in " << diff << " ms" << std::endl;
	}
}

void RaySpaceTree::constructAdaptive(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
									std::vector<std::pair<int,int>>& samples, std::set<int>& tris, bool print, std::vector<Ray> splitters) {
	auto start_time = std::chrono::high_resolution_clock::now();
	if (print) std::cout << "level: " << level << std::endl;
	bool checkmore = false;
	if (tris.size() == 1) {
		for (auto s : samples)
			if (s.second < 0) checkmore = true;
	}

	if (tris.size() > 1 || checkmore) {
		if (print) std::cout << "Node with " << samples.size() << " samples" << std::endl;
		if (print) std::cout << "Finding best splitter..." << std::endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		Ray splitter;
		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
		std::set<int> trisL = std::set<int>();
		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
		std::set<int> trisR = std::set<int>();
		if (getBestSplitter(splitter, res, samplesL, trisL, samplesR, trisR, cams, 10, samples, tris, print, splitters)) {

			auto end_time = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
			if (print) std::cout << "Found splitter in " << diff << " ms" << std::endl;

			node->splitter = splitter;
			splitters.push_back(splitter);

			Node* leftnode = new Node(node->index * 2 + 1, level + 1);
			node->leftNode = leftnode;
			leftnode->parent = node;
			nodes.push_back(leftnode);
			nodeSamples nsL = { trisL, samplesL, node->leftNode, level + 1 };
			toProcess.push(nsL);
			//constructAdaptive(level + 1, node->leftNode, res, cams, samplesL, trisL, print);

			Node* rightnode = new Node(node->index * 2 + 2, level + 1);
			rightnode->parent = node;
			node->rightNode = rightnode;
			nodes.push_back(rightnode);
			nodeSamples nsR = { trisR, samplesR, node->rightNode, level + 1 };
			toProcess.push(nsR);
			//constructAdaptive(level + 1, node->rightNode, res, cams, samplesR, trisR, print);

			return;
		}
	}
	noLeaves++;
	if (print) std::cout << "Leaf with " << samples.size() << " samples and " << tris.size() << " triangles" << std::endl;
	node->leaf = true;
	if (tris.size() > 0) {
		node->primitiveSet = tris;
		for (auto sample : samples) {
			if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
		if (print) std::cout << "Filled leaf in " << diff << " ms" << std::endl;
	}
}

Ray RaySpaceTree::getRay(std::vector<Camera*>& cams, int raynr, glm::ivec2 res) {
	int camNr= raynr / (res.x*res.y);

	int rayInCam = raynr % (res.x * res.y);
	int y = rayInCam / res.x;
	int x = rayInCam % res.x;
	const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
	return cams[camNr]->pixRayDirection(pixpos);
}

bool RaySpaceTree::intoLeftNode(Ray &splitter, std::vector<Camera*>& cams, int raynr, glm::ivec2 res) {
	Ray ray = getRay(cams, raynr, res);
	return splitter.side(ray);
}

Ray RaySpaceTree::getBestSplitter2(glm::ivec2 res, std::vector<Camera*>& cams, int tries,
									std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize) {
	int t = (int)time(NULL);
	srand(t);
	float bestTriFrac = 0;
	Ray splitter;

	//for (int i = 0; i < tris.size(); i++) {
		//for (int j = 0; j < 3; j++) {
	for (int i = 0; i < tries; i++) {
			std::cout << "Checking line no " << i << "... ";

			// choose the first edge of a random triangle in this node
			float r = (float)rand() / static_cast <float> (RAND_MAX);
			int rtri = int(r * (float)tris.size());
			float r2 = (float)rand() / static_cast <float> (RAND_MAX);
			int redge = int(r2 * 3.f);
			std::set<int>::iterator it = tris.begin();
			std::advance(it, rtri);
			int tri = *it;
			Ray testRay = Ray(model->vertices[3 * tri + redge].pos, model->vertices[3 * tri + (redge + 1) % 3].pos);
			//Ray testRay = Ray(vertices[3 * i + j].pos, vertices[3 * i + (j + 1) % 3].pos);


			std::set<int> trisLeft = std::set<int>();
			std::set<int> trisRight = std::set<int>();
			int zeroL = 0; int zeroR = 0;

			for (auto sample : samples) {
				bool left = intoLeftNode(testRay, cams, sample.first, res);
				int tri = sample.second;
				if (left) {
					if (tri >= 0) {
						trisLeft.insert(tri);

					}
					else zeroL++;
				}
				else {
					if (tri >= 0) {
						trisRight.insert(tri);
					}
					else zeroR++;
				}
			}

			float frac = 0;
			if (trisLeft.size() && trisRight.size()) frac = 1.f*trisLeft.size() / (1.f*trisRight.size());
			if (frac > 1.f) frac = 1.f / frac;
			std::cout << " for a tri frac of " << frac;

			if (frac > bestTriFrac) {
				std::cout << " which is better!" << std::endl;
				bestTriFrac = frac;
				splitter = testRay;
			}
			else {
				std::cout << std::endl;
			}
	//	}
	}
	std::cout << "New value: " << bestTriFrac << std::endl;
	//if (bestTriFrac == 0) getBestSplitter(splitter, res, cams, tries, samples, tris, totSampSize);
	return splitter;
}

bool RaySpaceTree::getBestSplitter(Ray& splitter, glm::ivec2 res, std::vector<std::pair<int, int>>& samplesLx, std::set<int>& trisLx, 
									std::vector<std::pair<int, int>>& samplesRx, std::set<int>& trisRx, std::vector<Camera*>& cams, int tries,
									std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray>& splitters) {
	int t = (int)time(NULL);
	srand(t);

	float bestSAH = tris.size() * 1.f;
	if (print) std::cout << "Parent SAH = " << bestSAH << std::endl;
	float parentSAH = bestSAH;

	std::vector<Ray> testRays = std::vector<Ray>();
	std::vector<int> triIndices = std::vector<int>();
	std::vector<int> triIndices2 = std::vector<int>();
	
	if (tris.size() <= 3) {
		for (int tri : tris) {
			for (int j = 0; j < 3; j++) {
				int triIndex = tri;
				int fromIndex = model->indices[3 * tri + j];
				int toIndex = model->indices[3 * tri + (j + 1) % 3];
				Ray testRay = Ray(model->vertices2[fromIndex], model->vertices2[toIndex]);
	
				//int triIndex2 = -5;
				//if (model->triPerVertex[fromIndex].size() > 1 && model->triPerVertex[toIndex].size() > 1) {
				//	for (auto trix : model->triPerVertex[fromIndex]) {
				//		if (triIndex2 >= 0) break;
				//		for (auto triy : model->triPerVertex[fromIndex]) {
				//			if (trix != triIndex && trix == triy) {
				//				triIndex2 = trix;
				//				break;
				//			}
				//		}
				//	}
				//}
				testRays.push_back(testRay);
				//triIndices.push_back(triIndex);
				//.push_back(triIndex2);
			}
		}
	}
	else {
		int i = 0;
		int upperlimit = 2 * tries;
		while (i < tries && i < upperlimit) {
			float r = (float)rand() / static_cast <float> (RAND_MAX);
			int rtri = int(r * (float)tris.size());
			float r2 = (float)rand() / static_cast <float> (RAND_MAX);
			int redge = int(r2 * 3.f);
			std::set<int>::iterator it = tris.begin();
			std::advance(it, rtri);
			int triIndex = *it;
			int fromIndex = model->indices[3 * triIndex + redge];
			int toIndex = model->indices[3 * triIndex + (redge + 1) % 3];
			Ray testRay = Ray(model->vertices2[fromIndex], model->vertices2[toIndex]);
			bool dupli = false;
			for (auto r : splitters) if (r.equal(testRay, 1E-6)) { dupli = true; break; }
			if (dupli) continue;
			i++;
			testRays.push_back(testRay);
	
			//int triIndex2 = -5;
			//if (model->triPerVertex[fromIndex].size() > 1 && model->triPerVertex[toIndex].size() > 1) {
			//	for (auto trix : model->triPerVertex[fromIndex]) {
			//		if (triIndex2 >= 0) break;
			//		for (auto triy : model->triPerVertex[fromIndex]) {
			//			if (trix != triIndex && trix == triy) {
			//				triIndex2 = trix;
			//				break;
			//			}
			//		}
			//	}
			//}	
			//triIndices.push_back(triIndex);
			//triIndices2.push_back(triIndex2);
		}
		if (testRays.size() == 0) std::cout << "Could not find any non-duplicate candidates" << std::endl;
	}

	// GET NODE AND THEN TRACE BACK WHICH SPLITTERS ARE ALREADY IN USE.
	// STORE IN SET OR IN VECTOR?
	// THEN ADD EVERYTIME SOMETHING NEW IS TRIED.
	// AND CHECK IF NEW OPTION IT'S ALREADY USED.

	// What if one of the samples incorrectly hits empty space? -- with the additional raytracing that should
	// now be the only check that needs to be done.
	for (int i = 0; i < testRays.size(); i++) {
		//std::vector<int> checktri = { triIndices[i], triIndices2[i] };
		//std::vector<std::vector<int>> checkIndTris = std::vector<std::vector<int>>(2*checktri.size());

		std::vector<std::pair<int, int>> samplesL = std::vector<std::pair<int, int>>();
		std::set<int> trisL = std::set<int>();
		std::vector<std::pair<int, int>> samplesR = std::vector<std::pair<int, int>>();
		std::set<int> trisR = std::set<int>();

		int pixL = 0;
		int pixR = 0;
		for (auto sample : samples) {
			bool left = intoLeftNode(testRays[i], cams, sample.first, res);
			int tri = sample.second;

			//for (int q = 0; q < checktri.size(); q++) {
			//	if (checktri[q] == tri) {
			//		if (left) checkIndTris[2*q].push_back(pixL);
			//		else checkIndTris[2*q+1].push_back(pixR);;
			//	}
			//}
			if (left) {
				samplesL.push_back(sample);
				if (tri >= 0) trisL.insert(tri);
				pixL++;
			}
			else {
				samplesR.push_back(sample);
				if (tri >= 0) trisR.insert(tri);
				pixR++;
			}
		}

		int samL = samplesL.size();
		int samR = samplesR.size();
		//std::vector<int> toRemoveL;
		//std::vector<int> toRemoveR;

		//for (int q = 0; q < checktri.size(); q++) {
		//	int leftsize = checkIndTris[2 * q].size();
		//	int rightsize = checkIndTris[2 * q + 1].size();
		//	if (leftsize > 0 && rightsize > 0) {
		//		if (leftsize < rightsize) {
		//			//for (auto sp : checkIndTris[2 * q]) {
		//			//	samplesR.push_back(samplesL[sp]);
		//			//	samplesL[sp] = samplesL.back();
		//			//	//samplesL.pop_back();
		//			//}
		//			samL--;
		//			samR++;
		//			trisL.erase(checktri[q]);
		//		}
		//		else {
		//			//for (auto sp : checkIndTris[2 * q + 1]) {
		//				//samplesL.push_back(samplesR[sp]);
		//				//samplesR[sp] = samplesR.back();
		//				//samplesR.pop_back();
		//			//}
		//			samL++;
		//			samR--;
		//			trisR.erase(checktri[q]);
		//		}
		//	}
		//}

		float fracLeft = 1.f * (samL) / float(samples.size());
		float fracRight = 1.f * (samR) / float(samples.size());

		float sahleft = trisL.size() * fracLeft; // sah
		float sahright = trisR.size() * fracRight;
		if (print) std::cout << "with " << fracLeft << ", " << fracRight << " samples (L,R) and " << trisL.size() << ", " << trisR.size() << " triangles (L,R) for a SAH of " << sahleft + sahright;

		if (sahleft + sahright < bestSAH) {
			if (print) std::cout << " which is better!" << std::endl;
			bestSAH = sahleft + sahright;
			splitter = testRays[i];
			samplesLx = samplesL;
			samplesRx = samplesR;
			trisLx = trisL;
			trisRx = trisR;
		}
		else {
			if (print) std::cout << std::endl;
		}
	}
	if (print) std::cout << "Old value: " << tris.size() << ", New value: " << bestSAH << std::endl;
	if (tris.size() == bestSAH) return false;
	return true;
}

void RaySpaceTree::printTree() {
	printTree(rootNode, 0);
	std::cout << std::endl << std::endl;
}

void RaySpaceTree::printTree(Node* node, int level) {
	if (!node->leaf) {
		level++;
		std::cout << std::endl;
		for (int x = 0; x < level; x++) std::cout << "  ";
		std::cout << "L";
		printTree(node->leftNode, level);
		std::cout << std::endl;
		for (int x = 0; x < level; x++) std::cout << "  ";
		std::cout << "R";
		printTree(node->rightNode, level);
	}
	else {
		std::cout << " Leaf: tri = " << node->primitiveSet.size() << ", samples = " << node->primAndRayVector.size();
	}
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, bool putRay, bool putPrim) {
	putPrimitive(ray, primId, rootNode, putRay, putPrim);
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, Node* node, bool putRay, bool putPrim) {
	if (node->leaf) {
		if (putRay) node->insert(primId, ray);
		else if (putPrim) node->insert(primId);
		return;
	}
	if (node->splitter.side(ray)) putPrimitive(ray, primId, node->leftNode, putRay, putPrim);
	else putPrimitive(ray, primId, node->rightNode, putRay, putPrim);
}

Node* RaySpaceTree::descend(Ray& ray)
{
	return descend(ray, rootNode);
}

Node* RaySpaceTree::descend(Ray& ray, Node* node)
{
	if (node->leaf) return node;
	if (node->splitter.side(ray)) return descend(ray, node->leftNode);
	else return descend(ray, node->rightNode);
}

std::vector<int> RaySpaceTree::countDuplicates(int size, std::vector<int>& nodenr) {
	std::vector<int> duplicates = std::vector<int>(size);
	countDuplicates(rootNode, duplicates, nodenr);
	return duplicates;
}

void RaySpaceTree::countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr) {
	if (!node->leaf) {
		countDuplicates(node->leftNode, duplicates, nodenr);
		countDuplicates(node->rightNode, duplicates, nodenr);
	}
	else {
		for (int x : node->primitiveSet) {
			duplicates[x] += 1;
			nodenr[node->index] += 1;
		}
	}
}

int RaySpaceTree::getNumberOfTriInleaf(int leafnum) {
	Node* n = getLeafFromNum(leafnum);
	return n->primitiveSet.size();
}

void RaySpaceTree::getSplittingLinesInLeaf(int leafnum, std::vector<Ray>& lines) {
	Node* leaf = getLeafFromNum(leafnum);
	getSplittingLinesInLeaf(leaf, lines);
};

void RaySpaceTree::getSplittingLinesInLeaf(Node* n, std::vector<Ray>& lines) {
	if (n == rootNode) return;
	else {
		lines.push_back(n->parent->splitter);
		return getSplittingLinesInLeaf(n->parent, lines);
	}
};

std::vector<Ray> RaySpaceTree::getAllSplittingLines() {
	std::vector<Ray> lines = std::vector<Ray>();
	for (Node* n : nodes) {
		if (n->leaf) continue;
		lines.push_back(n->splitter);
	}
	//getSplittingLines(rootNode, lines);
	return lines;
}

void RaySpaceTree::getSplittingLines(Node* node, std::vector<Ray>& lines) {
	if (!node->leaf) {
		lines.push_back(node->splitter);
		getSplittingLines(node->leftNode, lines);
		getSplittingLines(node->rightNode, lines);
	}
}

std::vector<glm::vec3> RaySpaceTree::getSplittingLinesInGeo(GeoObject* object) {

	std::vector<Ray> lines = getAllSplittingLines();
	std::vector<glm::vec3> splitters;
	glm::vec3 start, end;

	for (Ray line : lines) {
		object->intersect(line, start, end);
		splitters.push_back(start);
		splitters.push_back(end);
	}
	return splitters;
}

void RaySpaceTree::getSplittingLinesInLeafWithSide(Node* n, std::vector<Ray>& lines, std::vector<bool>& sides) {
	if (n == rootNode) return;
	else {
		if (n->parent->leftNode == n) sides.push_back(true);
		else sides.push_back(false);
		lines.push_back(n->parent->splitter);
		return getSplittingLinesInLeafWithSide(n->parent, lines, sides);
	}
} 

void RaySpaceTree::filterSplittingLines(Node* leaf, std::vector<Ray>& splitlines, std::vector<bool>& sides, 
										std::vector<Ray>& filteredLines, std::vector<bool>& filteredSides) {
	std::vector<int> foundExtremalStabbing;
	std::vector<std::vector<int>>& splitCombi4 = Combinations::combi4(splitlines.size());
	//for (std::vector<int>& combi4 : splitCombi4) {
	for (int i = 0; i < splitCombi4.size(); i++) {
		std::vector<Ray> lines4;
		for (int c : splitCombi4[i]) lines4.push_back(splitlines[c]);
		std::vector<Ray> intersectLines = LineThroughFour::find(lines4, model);
		for (Ray& r : intersectLines) {
			if (checkRayInLeaf(leaf, r, lines4, 4, false)) {
				foundExtremalStabbing.push_back(i);
				break;
			}
		}
	}

	for (int i = 0; i < splitlines.size(); i++) {
		bool found = false;
		for (int esl : foundExtremalStabbing) {
			if (found) break;
			for (int line : splitCombi4[esl]) {
				if (line == i) {
					filteredLines.push_back(splitlines[i]);
					if (i < sides.size()) filteredSides.push_back(sides[i]);
					found = true;
					break;
				}
			}
		}
	}
}

//bool RaySpaceTree::onCorrectSide(std::vector<Ray>& lines, std::vector<int>& sides, Ray& r) {
//	for (int i = 0; i < lines.size(); i++) {
//		if (lines[i].sideVal(r) > -1E-8 && sides[i])
//			return false;
//	}
//	return true;
//}

int RaySpaceTree::numOfLeaf(int index) {
	int count = 0;
	for (Node* n: nodes) {
		if (n->leaf) {
			count++;
			if (n->index == index) return count;
		}
	}
	return 0;
}

Node* RaySpaceTree::getLeafFromNum(int leafNum) {
	int leaf = 0;
	for (Node* n : nodes) {
		if (n->leaf) {
			if (leaf == leafNum) {
				return n;
			}
			leaf++;
		}
	}
	return NULL;
}

std::vector<Ray> RaySpaceTree::getViewingLinesInLeaf(int leafNum)
{
	Node* node = getLeafFromNum(leafNum);
	return getViewingLinesInLeaf(node);
}

std::vector<Ray> RaySpaceTree::getViewingLinesInLeaf(Node *node)
{
	std::vector<Ray> rays;
	for (Sample& p : node->primAndRayVector) rays.push_back(p.ray);
	return rays;
}

// add ignore list
bool RaySpaceTree::checkLineInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool& inPlane, bool print) {
	//check orientation
	float orient = glm::dot((glm::vec3)line.direction, model->normalPerTri[prim]);
	if (orient > 0) {
		if (print) std::cout << "Ray not in Prim (wrong orientation): " << orient << std::endl << std::endl;
		return false;
	}

	for (auto& er : edgeRays) {
		bool equal = false;
		for (auto& igray : lines4) {
			if (igray.equal(er, 1E-8)) { //-8
				equal = true;
				break;
			}
		}
		if (!equal && er.sideVal(line) < -1E-8) { //-8
			if (print) std::cout << "Ray not in Prim (wrong side of edge): " << er.sideVal(line) << std::endl << std::endl;
			return false;
		}
	}

	if (line.inPlane(edgeRays[0], edgeRays[1], 1E-10)) {
		if (print) std::cout << "Ray in plane: " << " ";
		for (int i = 0; i < 3; i++) {
			glm::dvec3 vert = model->vertices[3 * prim + i].pos;
			if (line.throughVertex(vert, 1E-10)) {
				if (print) std::cout << std::endl;
				inPlane = true;
				return true; 
			}
		}
		if (print) std::cout << std::endl;
		return false;
	}
	return true;
}

bool RaySpaceTree::checkRayInLeaf(Node* node, const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print) {

	if (node == rootNode) return true;
	Node* parent = node->parent;
	bool ign = false;
	for (int i = 0; i < rayIgnoresize; i++) {
		if (lines[i].equal(parent->splitter, 1E-8)) {
			ign = true; //-8
			break;
		}
	}
	double sidev = parent->splitter.sideVal(ray);// (was other way around!) .sideVal;
	if (abs(sidev) < 1E-8) ign = true;

	if (parent->leftNode == node) {
		if (sidev < 0 || ign) return checkRayInLeaf(parent, ray, lines, rayIgnoresize, print);
		if (print) std::cout << "Ray not in Leaf (on wrong side of splitting line, should be <0): " << sidev << std::endl << std::endl;
	}
	else {
		if (sidev > 0 || ign) return checkRayInLeaf(parent, ray, lines, rayIgnoresize, print);
		if (print) std::cout << "Ray not in Leaf: on wrong side of splitting line, should be >0 " << sidev << std::endl << std::endl;
	}
	return false;
}

bool RaySpaceTree::checkPrimVisibleForLine(Ray& ray, const int prim, std::vector<int>& ignore, std::vector<Ray>& lines,  bool inplane, bool print) {

	int embreePrim = -1;
	float embreeDepth = 0.f;
	float primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
	//model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
	//if (embreePrim == prim || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
	//if (inplane) {
	//	glm::vec3 v1 = model->vertices[3 * embreePrim].pos;
	//	glm::vec3 v2 = model->vertices[3 * embreePrim + 1].pos;
	//	glm::vec3 v3 = model->vertices[3 * embreePrim + 2].pos;
	//	std::vector<Ray> edgeRays = { Ray(v2, v1), Ray(v3, v2), Ray(v1, v3) };
	//	glm::dvec3 cross = glm::cross(glm::normalize(edgeRays[0].direction), glm::normalize(edgeRays[1].direction));
	//	double v = fabsf(glm::dot(ray.direction, cross));
	//	if (v < 1E-10) return true;
	//}

	if (ignore.size() > 0) {
		std::vector<float> uniquePrimDepth;

		// goes wrong right now
		for (int i = 0; i < ignore.size() / 2; i++) {
			float tryDepth = model->getIntersectionDepthForPrim(ignore[i * 2] * 3, ray);
			bool found = false;
			for (int j = 0; j < uniquePrimDepth.size(); j++) {
				if (fabsf(tryDepth - uniquePrimDepth[j]) < 1E-3) {
					found = true;
					break;
				}
			}
			if (!found) uniquePrimDepth.push_back(tryDepth);
		}

		int i = 0;
		while (i < uniquePrimDepth.size()) {
		//for (int i = 0; i < uniquePrimDepth.size(); i++) {
			bool found = false;
			model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
			for (int j = 0; j < ignore.size() / 2; j++) {
				// only checking first ignore primitive for now
				float primdepth = model->getIntersectionDepthForPrim(ignore[j * 2] * 3, ray);
				// check ray direction vs two triangles to see if it is still silhouette
				if (ignore[j*2+1] >= 0)
					if ((glm::dot(model->normalPerTri[ignore[j * 2]], (glm::vec3)ray.direction) < 0) ==
					(glm::dot(model->normalPerTri[ignore[j * 2 + 1]], (glm::vec3)ray.direction) < 0)) continue;
				
				// Checking primary primdepth here does not result in 'correct' Extremal Stabbing Line
				primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
				if (fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
				if (fabsf(embreeDepth - primdepth) < 1E-3) {
					glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
					ray = Ray(neworig + ray.direction, neworig);
					found = true;
					break;
				}
			}
			if (!found) {
				if (lines.size() > 0) {
					for (Ray& r : lines) {
						glm::dvec3 intersectionTs = (ray.origin + ray.direction * (double)embreeDepth - r.origin) / r.direction;
						if (fabsf(intersectionTs.x - intersectionTs.y) < 1E-10) {
							glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
							ray = Ray(neworig + ray.direction, neworig);
							i--;
							continue;
						}
					}
				}
				return false;
			}
			i++;
		}
	}

	embreePrim = -1;
	embreeDepth = 0.f;
	primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
	model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
	if (embreePrim >= 0) {
		if (embreePrim == prim || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
		else if (inplane) {
			glm::vec3 v1 = model->vertices[3 * embreePrim].pos;
			glm::vec3 v2 = model->vertices[3 * embreePrim + 1].pos;
			glm::vec3 v3 = model->vertices[3 * embreePrim + 2].pos;
			std::vector<Ray> edgeRays = { Ray(v2, v1), Ray(v3, v2), Ray(v1, v3) };
			glm::dvec3 cross = glm::cross(glm::normalize(edgeRays[0].direction), glm::normalize(edgeRays[1].direction));
			double v = fabsf(glm::dot(ray.direction, cross));
			if (v < 1E-10) return true;
		}
		else if (lines.size() > 0) {
			for (Ray& r : lines) {
				glm::dvec3 intersectionTs = (ray.origin + ray.direction * (double)embreeDepth - r.origin) / r.direction;
				if (fabsf(intersectionTs.x - intersectionTs.y) < 1E-10) {
					glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
					ray = Ray(neworig + ray.direction, neworig);
					primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
					model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
					if (embreePrim == prim || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
				}
			}
		}
	}
	return false;
};

bool RaySpaceTree::checkLineInBox(const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print) {
	return model->boundingCube.intersectSide(maindir, ray, lines, rayIgnoresize, print);
}

//bool RaySpaceTree::findExtremalStabbingForPrim(const int prim, std::vector<std::vector<int>>& combi, std::vector<Ray>& splitLines,
//												Ray& ray, Node* leaf, std::vector<Ray>& edgeRays, bool print, int nrOfsilhEdges,
//												std::vector<Ray>& silhLines, 
//												std::vector<int>& visibleTriIgnore, std::vector<Edge>& silhouetteEdges) {
//
//	std::vector<std::vector<int>> edges = { {2,0}, {0,1}, {1,2} };
//	int edgeCombi = 3;
//	int lineCombiSize = combi.size() > 0 ? combi[0].size() : 0;
//	if (lineCombiSize == 4) edgeCombi = 1;
//	int nrOfSplitLines = lineCombiSize - nrOfsilhEdges;
//
//	for (int e = 0; e < edgeCombi; e++) {
//		for (int c = 0; c < combi.size(); c++) {
//
//			std::vector<Ray> lines;
//			for (int i = 0; i < nrOfSplitLines; i++) lines.push_back(splitLines[combi[c][i]]);
//			for (int i = 0; i < nrOfsilhEdges; i++) lines.push_back(silhLines[combi[c][i]]);
//			for (int i = 0; i < 4 - lineCombiSize; i++) lines.push_back(edgeRays[edges[e][i]]);
//			
//			if (nrOfsilhEdges >= 2 && nrOfsilhEdges < 4) {
//
//			}
//
//			std::vector<int> triIgnore;
//			for (int i = 0; i < nrOfsilhEdges; i++) {
//				triIgnore.push_back(visibleTriIgnore[(combi[c][i]) * 2]);
//				triIgnore.push_back(visibleTriIgnore[(combi[c][i]) * 2 + 1]);
//			}
//			if (checkRaysThroughLines(prim, ray, leaf, nrOfSplitLines, edgeRays, print, triIgnore, lines, lineCombiSize)) return true;
//			
//		}
//	}
//	return false;
//}
bool RaySpaceTree::checkCombi(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::string combi_text, int combiNr,
								int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges,
								std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis,
								std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& silhLineCombis, std::vector<Edge>& silhouetteEdges,
								std::vector<Ray>& silhVertexLines, std::vector<std::vector<int>>& silhVertexCombis,
								std::vector<int>& silhouetteTris, std::vector<Ray>& triEdgeRays) {

	if (print) std::cout << combi_text << " combi's: " << combiNr << std::endl;

	std::vector<std::vector<int>> edges = { {2,0}, {0,1}, {1,2} };
	std::vector<int> vertexEdgeCheck;

	for (int g = 0; g < std::min(1, nrOfTriEdges*3); g++) {
		if (nrOfTriEdges == 2) vertexEdgeCheck.push_back(model->indices[prim * 3 + edges[g][0]]);
		if (nrOfTriEdges == 1) vertexEdgeCheck.push_back(model->indices[prim * 3 + edges[g][1]]);

		for (int h = 0; h < std::min(1, (int)silhVertexCombis.size()); h++) {
			if (nrOfVertices > 0) vertexEdgeCheck.push_back(silhVertexCombis[h][3]);

			for (int i = 0; i < std::min(1, (int)silhLineCombis.size()); i++) {
				bool cntn;
				for (int ix = 0; ix < nrOfsilhEdges; ix++) {
					cntn = true;
					if (nrOfsilhEdges > 0) {
						Edge e = silhouetteEdges[silhLineCombis[i][ix]];
						if (vertexEdgeCheck.size() == 1) {
							if (vertexEdgeCheck[0] == *e.vertices.begin() || vertexEdgeCheck[0] == *e.vertices.rbegin()) break;
							bool side;
							if (model->checkSilhouetteEdge2(model->verticesIndexed[vertexEdgeCheck[0]], e, true, glm::vec3(0), side) <= 0) break;
						}
					}
					if ((nrOfsilhEdges > 2 || (nrOfsilhEdges == 2 && ix == 0)) && vertexEdgeCheck.size() > 0) {
						std::vector<glm::vec3> n;
						std::vector<float> d;
						spaceSpannedByEdges(silhouetteEdges[silhLineCombis[i][ix]], silhouetteEdges[silhLineCombis[i][(ix + 1) % nrOfsilhEdges]], n, d);
						std::vector<glm::vec3> checkPoints;
						for (int vec : vertexEdgeCheck) checkPoints.push_back(model->verticesIndexed[vec]);
						if (!checkPointsInHalfSpaces(n, d, checkPoints)) break;
					}
					cntn = false;
				}
				if (cntn) continue;

				for (int j = 0; j < std::min(1, (int)splitLineCombis.size()); j++) {

					std::vector<Ray> lines;
					std::vector<int> tris;

					for (int k = 0; k < nrOfsplitLines; k++) lines.push_back(splitLines[splitLineCombis[j][k]]);
					if (nrOfVertices) {
						lines.push_back(silhouetteLines[silhLineCombis[h][0]]);
						lines.push_back(silhVertexLines[silhLineCombis[h][1]]);
						tris.push_back(silhouetteTris[2 * silhLineCombis[h][0]]);
						tris.push_back(silhouetteTris[2 * silhLineCombis[h][0] + 1]);
					}
					for (int k = 0; k < nrOfsilhEdges; k++) {
						lines.push_back(silhouetteLines[splitLineCombis[i][k]]);
						tris.push_back(silhouetteTris[2 * splitLineCombis[i][k]]);
						tris.push_back(silhouetteTris[2 * splitLineCombis[i][k] + 1]);
					}
					for (int k = 0; k < nrOfTriEdges; k++) lines.push_back(triEdgeRays[edges[k][i]]);

					if (checkRaysThroughLines(prim, ray, leaf, 0, triEdgeRays, printAll, tris, lines)) {
						if (print) std::cout << combi_text << std::endl;
						return true;
					}
				}
			}
		}
	}
}



bool RaySpaceTree::checkRayThroughVertices(int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
												bool print, std::vector<int>& triIgnore, std::vector<Ray>& lines, int rayIgnoresize) {
	for (int i = 0; i < 3; i++) {
		if (ray.throughVertex(model->vertices[prim * 3 + i].pos, 1E-7)) {
			bool found = false;
			Vertex v = model->vertices[prim * 3 + i];
			glm::vec3 normal = model->normalPerTri[prim];
			glm::vec3 from = v.pos + (v.center - v.pos) * 0.1f;
			glm::vec3 to0 = v.pos + (glm::vec3)ray.direction;
			Ray ray0 = Ray(to0, from);
			found = checkExtremalStabbingLine(prim, ray0, leaf, splitsize, edgeRays, print, triIgnore, lines, rayIgnoresize);
			if (!found) {
				glm::vec3 to1 = v.pos + (glm::vec3)ray.direction + 0.1f * normal;
				Ray ray1 = Ray(to1, from);
				found = checkExtremalStabbingLine(prim, ray1, leaf, splitsize, edgeRays, print, triIgnore, lines, rayIgnoresize);
			}
			if (!found) {
				glm::vec3 to2 = (glm::vec3)ray.direction - 0.1f * normal;
				Ray ray2 = Ray(to2, from);
				found = checkExtremalStabbingLine(prim, ray2, leaf, splitsize, edgeRays, print, triIgnore, lines, rayIgnoresize);
			}
			if (!found) return false;
		}
	}
}

bool RaySpaceTree::checkExtremalStabbingLine(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
												bool print, std::vector<int>& triIgnore, std::vector<Ray>& lines, int rayIgnoresize) {

	if (!alldir) if (!checkLineInBox(ray, lines, rayIgnoresize, print)) return false;
	if (!checkRayInLeaf(leaf, ray, lines, rayIgnoresize, print)) return false;

	// not completely correct, should test for specific splitting line that also goes through this vertex.
	ray.get3DfromPlucker();
	bool inPlane = false;
	if (!checkLineInPrim(edgeRays, ray, lines, prim, inPlane, print)) return false;

	// not for generic viewing direction yet, only the three 'minima'
	double ttest = glm::dot(((model->boundingBox.getBounds(0) - glm::vec3(0.1) - (glm::vec3)ray.origin) / (glm::vec3)ray.direction), glm::abs(maindir));
	glm::dvec3 neworig = glm::dvec3(ray.origin + ttest * ray.direction);
	Ray ray2 = Ray(neworig + ray.direction, neworig);
	if (!checkPrimVisibleForLine(ray2, prim, triIgnore, lines, inPlane, print)) return false;

	//if (!checkRayThroughVertices(prim, ray, leaf, splitsize, edgeRays, print, triIgnore, lines, rayIgnoresize)) return false;

	//if (print && triIgnore.size() == 2) {
	//	std::cout << "Silhouette Triangles = " << triIgnore[0] << " & " << triIgnore[1] << std::endl;
	//}

	//if (print) {
	//	std::cout << std::endl << " found ray ";
	//	ray.print();
	//	std::cout << std::endl;
	//}
	return true;
}


bool RaySpaceTree::checkRaysThroughLines(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
											bool print, std::vector<int>& visibleTriIgnore, std::vector<Ray>& lines, int rayIgnoresize) {
	std::vector<Ray> intersectLines = LineThroughFour::find(lines, model);
	for (int i = 0; i < intersectLines.size(); i++) {
		ray = intersectLines[i];
		if (checkRayAndReverse(prim, ray, leaf, 0, edgeRays, print, visibleTriIgnore)) return true;
	}
	return false;
}

bool RaySpaceTree::checkRayAndReverse(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
									bool print, std::vector<int>& visibleTriIgnore, std::vector<Ray>& lines, int rayIgnoresize) {
	if (checkExtremalStabbingLine(prim, ray, leaf, 0, edgeRays, print, visibleTriIgnore)) return true;
	ray.inverseDir();
	if (checkExtremalStabbingLine(prim, ray, leaf, 0, edgeRays, print, visibleTriIgnore)) return true;
	return false;
}

bool RaySpaceTree::check1Prim(const int prim, Ray& ray, Node* leaf, bool print, int edgeSelection) {
	std::vector<Ray> splitLines;
	std::vector<bool> sideLines;

	getSplittingLinesInLeafWithSide(leaf, splitLines, sideLines);
	if (!alldir) {
		std::vector<Ray> boxSides = model->boundingCube.getCubeSideLines(maindir);
		for (auto r : boxSides) splitLines.push_back(r);
	}
	int size = splitLines.size();// +boxSides.size();
	std::vector<std::vector<int>> combi2 = Combinations::combi2(size);
	std::vector<std::vector<int>> combi3 = Combinations::combi3(size);
	std::vector<std::vector<int>> combi4 = Combinations::combi4(size);

	return checkPrim(prim, combi2, combi3, combi4, splitLines, sideLines, ray, leaf, print, edgeSelection);
}


//bool RaySpaceTree::checkSSVt(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<std::vector<int>>& splitCombi,
//							std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "SSV(t) combi's: " << splitCombi.size() * 3 << std::endl;
//	if (findExtremalStabbingForPrim(prim, splitCombi, splitLines, ray, leaf, edgeRays, printAll)) {
//		if (print) std::cout << "SSV(t)" << std::endl;
//		return true;
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSSST(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<std::vector<int>>& splitCombi, 
//							std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "SSST combi's: " << splitCombi.size() * 3 << std::endl;
//	if (findExtremalStabbingForPrim(prim, splitCombi, splitLines, ray, leaf, edgeRays, printAll)) {
//		if (print) std::cout << "SSST" << std::endl;
//		return true;
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSSSS(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<std::vector<int>>& splitCombi, 
//							std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "SSSS combi's: " << splitCombi.size() << std::endl;
//	if (findExtremalStabbingForPrim(prim, splitCombi, splitLines, ray, leaf, edgeRays, printAll)) {
//		if (print) std::cout << "SSSS" << std::endl;
//		return true;
//	}
//	return false;
//}
//
//
//
bool RaySpaceTree::checkVeVt(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges, 
							std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
	if (print) std::cout << "V(e)V(t) combi's: " << silhouetteEdges.size() * 2 * 3 << std::endl;
	// there are doubles in here! filter out to get unique vertices? a set perhaps??
	for (int i = 0; i < silhouetteEdges.size(); i++) {
		for (int v : silhouetteEdges[i].vertices) {
			std::vector<int> tris = { silhouetteTris[i * 2], silhouetteTris[i * 2 + 1] };
			for (int k = 0; k < 3; k++) {
				ray = Ray(model->verticesIndexed[v], model->vertices[prim * 3 + k].pos);
				if (checkRayAndReverse(prim, ray, leaf, 0, edgeRays, printAll, tris)) {
					if (print) std::cout << "V(e)V(t)" << std::endl;
					return true;
				}
			}
		}
	}
	return false;
}
//
//bool RaySpaceTree::checkSVeT(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitlines, 
//							std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, 
//							std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "SV(e)T combi's: " << splitlines.size() * silhouetteLines.size() * 2 * 3 << std::endl;
//	for (int i = 0; i < silhouetteLines.size(); i++) {
//		Edge e = silhouetteEdges[i];
//		Ray r1v = silhouetteLines[i];
//		std::vector<int> tris = { silhouetteTris[2 * i], silhouetteTris[2 * i + 1] };
//		for (int v : e.vertices) {
//			Ray r2v = Ray(model->vertices2[v] + model->normalPerTri[e.triangles[0]], model->vertices2[v]);
//			for (int s = 0; s < splitlines.size(); s++) for (Ray& t : edgeRays) {
//				if (checkRaysThroughLines(prim, ray, leaf, 0, edgeRays, printAll, tris, std::vector<Ray> { splitlines[s], r1v, r2v, t })) {
//					if (print) std::cout << "SV(e)T" << std::endl;
//					return true;
//				}
//			}
//		}
//	}
//	return false;
//}

//bool RaySpaceTree::checkVeET(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges, 
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "V(e)ET combi's: " << silhouetteLines.size() * silhouetteLines.size() * 2 * 3 << std::endl;
//	std::vector<int> tris(4);
//	for (int t = 0; t < 3; t++) {
//		Edge triEdge = { std::set<int>{model->vertices[prim * 3 + t].pos, model->vertices[prim * 3 + t % 3].pos} };
//		for (int i = 0; i < silhouetteLines.size(); i++) {
//			tris[0] = silhouetteTris[i * 2];
//			tris[1] = silhouetteTris[i * 2 + 1];
//			std::vector<glm::vec3> n;
//			std::vector<float> d;
//			spaceSpannedByEdges(silhouetteEdges[i], triEdge, n, d);
//
//			for (int j = i+1; j < silhouetteLines.size(); j++) {
//				Edge e = silhouetteEdges[j];
//				tris[2] = silhouetteTris[j * 2];
//				tris[3] = silhouetteTris[j * 2 + 1];
//				for (int v : e.vertices) {
//					if (!checkPointsInHalfSpaces(n, d, std::vector<glm::vec3>{model->vertices2[v]})) continue;
//					bool side;
//					if (model->checkSilhouetteEdge2(model->vertices2[v], silhouetteEdges[i], true, glm::vec3(0), side) <= 0) continue;
//
//					Ray r2v = Ray(model->vertices2[v] + model->normalPerTri[e.triangles[0]], model->vertices2[v]);
//
//					if (checkRaysThroughLines(prim, ray, leaf, 0, edgeRays, printAll, tris, std::vector<Ray> { silhouetteLines[i], silhouetteLines[j], r2v, edgeRays[t] })) {
//						if (print) std::cout << "V(e)ET" << std::endl;
//						return true;
//					}
//				}
//			}
//		}
//	}
//	return false;
//}

//bool RaySpaceTree::checkMutualEdgeVisibility() {
//	bool check12;
//	bool side1;
//	int silhEdgeCheck1 = model->checkSilhouetteEdge2(model->vertices2[*e1.vertices.begin()], e2, false, maindir, side1);
//	if (silhEdgeCheck1 == 0);
//	else if (silhEdgeCheck1 == 1) check12 = true;
//	else return false;
//	bool side2;
//	int silhEdgeCheck2 = model->checkSilhouetteEdge2(model->vertices2[*e1.vertices.rbegin()], e2, false, maindir, side2);
//	if (silhEdgeCheck1 == 1 || (silhEdgeCheck2 == 0 && side1 != side2)) check12 = true;
//	else return false;
//}


//bool RaySpaceTree::checkSVeE(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines,
//							std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, 
//							std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "SV(e)E combi's: " << splitLines.size() * silhouetteLines.size() * silhouetteLines.size() * 2 << std::endl;
//	std::vector<int> tris(4);
//	for (int i = 0; i < silhouetteLines.size(); i++) {
//		Edge e1 = silhouetteEdges[i];
//		tris[0] = silhouetteTris[i * 2];
//		tris[1] = silhouetteTris[i * 2 + 1];
//		for (int j = i + 1; j < silhouetteLines.size(); j++) {
//			Edge e2 = silhouetteEdges[j];
//			tris[2] = silhouetteTris[j * 2];
//			tris[3] = silhouetteTris[j * 2 + 1];
//			Ray r1v = silhouetteLines[j];
//			for (int v : e2.vertices) {
//				bool side;
//				if (model->checkSilhouetteEdge2(model->vertices2[v], e2, true, glm::vec3(0), side) <= 0) continue;
//				Ray r2v = Ray(model->vertices2[v] + model->normalPerTri[e2.triangles[0]], model->vertices2[v]);
//				for (int s = 0; s < splitLines.size(); s++) {
//					if (checkRaysThroughLines(prim, ray, leaf, 0, edgeRays, printAll, tris, std::vector<Ray> { splitLines[s], r1v, r2v, silhouetteLines[i]})) {
//						if (print) std::cout << "SV(e)E" << std::endl;
//						return true;
//					}
//				}
//			}
//		}
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSSVe(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<std::vector<int>>& splitCombi, 
//							std::vector<Ray>& splitLines, std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines,
//							std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//	if (print) std::cout << "SSV(e) combi's: " << splitLines.size() * splitLines.size() * silhouetteLines.size() * 2 << std::endl;
//	for (int i = 0; i < silhouetteLines.size(); i++) {
//		Edge e = silhouetteEdges[i];
//		std::vector<Ray> intersectLines;
//		Ray r1v = silhouetteLines[i];
//
//		std::vector<int> tris = { silhouetteTris[2 * i], silhouetteTris[2 * i + 1] };
//
//		for (int v : e.vertices) {
//			Ray r2v = Ray(model->vertices2[v] + model->normalPerTri[e.triangles[0]], model->vertices2[v]);
//			for (std::vector<int>& c: splitCombi){
//				if (checkRaysThroughLines(prim, ray, leaf, 0, edgeRays, printAll, tris, std::vector<Ray> { splitLines[c[0]], splitLines[c[1]], r1v, r2v})) {
//					if (print) std::cout << "SSV(e)" << std::endl;
//					return true;
//				}
//			}
//		}
//	}
//	return false;
//}




//bool RaySpaceTree::checkSEVt(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines, 
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//	std::vector<std::vector<int>> combi = Combinations::combi11(splitLines.size(), silhouetteLines.size());
//	if (print) std::cout << "SEV(t) combi's: " << combi.size() * 3 << std::endl;
//	if (findExtremalStabbingForPrim(prim, combi, splitLines, ray, leaf, edgeRays, printAll, 1, silhouetteTris, silhouetteLines)) {
//		if (print) std::cout << "SEV(t)" << std::endl;
//		return true;
//	}
//	return false;
//}

bool RaySpaceTree::checkVeVe(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges, 
							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
	//std::vector<std::vector<int>> combi = Combinations::combi2(silhouettesize);

	if (print) std::cout << "V(e)V(e) combi's: " << silhouetteLines.size() * 2 * silhouetteLines.size() * 2 << std::endl;
	std::vector<int> tris(4);
	for (int i = 0; i < silhouetteEdges.size(); i++) {
		Edge e1 = silhouetteEdges[i];
		tris[0] = silhouetteTris[i * 2];
		tris[1] = silhouetteTris[i * 2 + 1];
		for (int v1 : e1.vertices) {
			for (int j = i+1; j < silhouetteEdges.size(); j++) {
				Edge e2 = silhouetteEdges[j];
				tris[2] = silhouetteTris[j * 2];
				tris[3] = silhouetteTris[j * 2 + 1];

				for (int v2 : e2.vertices) {
					ray = Ray(model->verticesIndexed[v1], model->verticesIndexed[v2]);
					if (checkRayAndReverse(prim, ray, leaf, 0, edgeRays, printAll, tris)) {
						if (print) std::cout << "V(e)V(e)" << std::endl;
						return true;
					}
				}
			}
		}
	}
	return false;
}

void RaySpaceTree::spaceSpannedByEdges(Edge& e1, Edge& e2, std::vector<glm::vec3>& n, std::vector<float>& c, std::vector<glm::vec3>& triPlanes) {
	std::vector<glm::vec3> e1v = { model->vertices2[*e1.vertices.begin()], model->vertices2[*e1.vertices.rbegin()] };
	std::vector<glm::vec3> e2v = { model->vertices2[*e2.vertices.begin()], model->vertices2[*e2.vertices.rbegin()] };

	glm::vec3 e1v01 = e1v[0] - e1v[1];
	glm::vec3 e2v01 = e2v[0] - e2v[1];
	glm::vec3 e1v0_e2v0 = e1v[0] - e2v[0];
	glm::vec3 e1v1_e2v1 = e1v[1] - e2v[1];

	// 4 planes e1v0_e2v0-e1, e1v1_e2v1-e1, e1v0_e2v0-e2, e1v1_e2v1-e2
	n.resize(4);
	c.resize(4);
	triPlanes = {e1v[0], e1v[1], e2v[0], e1v[0], e1v[1], e2v[1], e2v[0], e2v[1], e1v[0], e2v[0], e2v[1], e1v[1] };

	addPlaneInVector(e1v01, e1v0_e2v0, e1v[0], e2v[1], n, c, 0);
	addPlaneInVector(e1v01, e1v1_e2v1, e1v[0], e2v[0], n, c, 1);
	addPlaneInVector(e2v01, e1v0_e2v0, e2v[0], e1v[1], n, c, 2);
	addPlaneInVector(e2v01, e1v1_e2v1, e2v[0], e1v[0], n, c, 3);
}

bool RaySpaceTree::checkPointsInHalfSpaces(std::vector<glm::vec3>& n, std::vector<float>& c, std::vector<glm::vec3>& points, float err) {
	// check if one of the two endpoints of the edge 3 is in the volume spanned by edge 1 and 2
	// 3 configuration, inside the box, side 1 and side 2
	// lots of double checks
	bool inVolume = true;
	for (int v = 0; v < 3; v++) {
		for (int i = 0; i < 4; i++) {
			if (v * 2 - i > 0) n[i] *= -1; // for v=1 only the first two flip, for v=2 everything flips
			int count = 0;
			for (glm::vec3& p : points) if (glm::dot(n[i], p) - c[i] > err) count++;
			if (count == points.size()) return false;
		}
		if (inVolume) return true;
	}
	return false;
}


//bool RaySpaceTree::checkVeEE(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges, 
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris,
//							std::vector<Ray>& edgeRays, std::vector<std::vector<int>>& combi2Edges) {
//	if (print) std::cout << "V(e)EE combi's: " << silhouetteLines.size() * 2 * combi2Edges.size() << std::endl;
//	for (std::vector<int>& c : combi2Edges) {
//		for (int i = 0; i < silhouetteEdges.size(); i++) {
//			if (i == c[0] || i == c[1]) continue;
//			Edge e = silhouetteEdges[i];
//			Edge e1 = silhouetteEdges[c[0]];
//			Edge e2 = silhouetteEdges[c[1]];
//			Ray r1v = silhouetteLines[i];
//
//			std::vector<glm::vec3> n;
//			std::vector<float> d;
//			spaceSpannedByEdges(e1, e2, n, d);
//			std::vector<int> tris = e.triangles;
//
//			for (int v : e.vertices) {
//				if (!checkPointsInHalfSpaces(n, d, std::vector<glm::vec3> {model->vertices2[v]})) continue;
//
//				// check vertex with both edges
//				Ray r2v = Ray(model->vertices2[v] + model->normalPerTri[e.triangles[0]], model->vertices2[v]);
//				if (checkRaysThroughLines(prim, ray, leaf, 0, edgeRays, printAll, tris, std::vector<Ray> { r1v, r2v, silhouetteLines[c[0]], silhouetteLines[c[1]]})) {
//					if (print) std::cout << "V(e)EE" << std::endl;
//					return true;
//				}
//			}
//		}
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSEET(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines,
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays,
//							std::vector<std::vector<int>>& allowedCombi) {
//	if (silhouetteLines.size() >= 2 && allowedCombi.size() > 0) {
//		std::vector<std::vector<int>> combi = Combinations::combiAdd1(splitLines.size(), allowedCombi);
//		if (print) std::cout << "SEET combi's: " << combi.size() * 3 << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi, splitLines, ray, leaf, edgeRays, printAll, 2, silhouetteTris)) {
//			if (print) std::cout << "SEET" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSSET(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines,
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//	if (silhouetteLines.size() >= 1) {
//		std::vector < std::vector<int>> combi = Combinations::combi21(splitLines.size(), silhouetteLines.size());
//		if (print) std::cout << "SSET combi's: " << combi.size() * 3 << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi, splitLines, ray, leaf, edgeRays, printAll, 1, silhouetteTris, silhouetteLines)) {
//			if (print) std::cout << "SSET" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSSSE(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines, 
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//	if (silhouetteLines.size() >= 1) {
//		std::vector < std::vector<int>> combi = Combinations::combi31(splitLines.size(), silhouetteLines.size());
//		if (print) std::cout << "SSSE combi's: " << combi.size() * 3 << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi, splitLines, ray, leaf, edgeRays, printAll, 1, silhouetteTris, silhouetteLines)) {
//			if (print) std::cout << "SSSE" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkSSEE(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines,
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays, 
//							std::vector<std::vector<int>>& allowedCombi) {
//	if (silhouetteLines.size() >= 2 && allowedCombi.size() > 0) {
//		// SS should intersect EEvolume
//		std::vector<std::vector<int>> combi = Combinations::combiAdd2(splitLines.size(), allowedCombi);
//		if (print) std::cout << "SSEE combi's: " << combi.size() << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi, splitLines, ray, leaf, edgeRays, printAll, 2, silhouetteTris, silhouetteLines)) {
//			if (print) std::cout << "SSEE" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}

void addPlaneInVector(glm::vec3 line1, glm::vec3 line2, glm::vec3 pointOnPlane, glm::vec3 pointNegativeSide, std::vector<glm::vec3>& n, std::vector<float>& c, int num) {

	n[num] = glm::cross(line1, line2);
	c[num] = glm::dot(n[num], pointOnPlane);
	if (glm::dot(n[num], pointNegativeSide) - c[num] > 0) n[num] *= -1;
}

//bool RaySpaceTree::checkEEET(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges, 
//							std::vector<Ray>& silhouetteLines,  std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays, 
//							std::vector<std::vector<int>>& allowedCombi, std::vector<std::vector<int>>& combi3Edges) {
//
void RaySpaceTree::getEEECombis(std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines,
	std::vector<std::vector<int>>& combi2Edges, std::vector<std::vector<int>>& combi3Edges) {
	if (silhouetteLines.size() >= 3 && combi2Edges.size() > 0) {
		std::vector<std::vector<int>> combi = Combinations::combiAddSelective(silhouetteLines.size(), combi2Edges);
		for (std::vector<int>& c : combi) {

			std::vector<glm::vec3> n;
			std::vector<float> d;
			std::vector<glm::vec3> tri;
			Edge e1 = silhouetteEdges[c[0]];
			Edge e2 = silhouetteEdges[c[1]];
			spaceSpannedByEdges(e1, e2, n, d, tri);
			std::vector<glm::vec3> e3 = { model->verticesIndexed[*silhouetteEdges[c[2]].vertices.begin()],
										  model->verticesIndexed[*silhouetteEdges[c[2]].vertices.rbegin()] };
			if (checkPointsInHalfSpaces(n, d, e3)) {
				combi3Edges.push_back(c);
				continue;
			}

			// find the intersection point of edge 3 with the planes spanned by edge 1 and 2
			for (int i = 0; i < 4; i++) {
				float t = (-glm::dot(e3[0], n[i]) + c[i]) / glm::dot(e3[1] - e3[0], n[i]);
				// if point does not lie between e3v0 and e3v1 discard
				if (t < 0 || t > 1) continue;

				// if point does not lie on triangle wedge (between edges) discard
				glm::vec3 is = e3[0] + (e3[1] - e3[0]) * t;
				bool dot1 = glm::dot(tri[i * 3 + 2] - tri[i * 3], is) > 0;
				bool dot2 = glm::dot(tri[i * 3 + 2] - tri[i * 3 + 1], is) > 0;
				if (dot1 != dot2) continue;

				// for now just check if it is in one of 3 parts
				if (checkPointsInHalfSpaces(n, d, std::vector<glm::vec3>{is}, 1E-8)) {
					combi3Edges.push_back(c);
					continue;
				}
			}
		}
	}
}
//		 
//		if (print) std::cout << "EEET combi's: " << combi3Edges.size() * 3 << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi3Edges, silhouetteLines, ray, leaf, edgeRays, printAll, 3, silhouetteTris)) {
//			if (print) std::cout << "EEET" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}
//
//bool RaySpaceTree::checkEEEE(const int prim, Ray & ray, Node * leaf, bool print, bool printAll, std::vector<Ray>&silhouetteLines, 
//							std::vector<int>&silhouetteTris, std::vector<Ray>&edgeRays, std::vector<std::vector<int>>& allowedCombi) {
//	if (silhouetteLines.size() >= 4 && allowedCombi.size() > 0) {
//		std::vector<std::vector<int>> combi = Combinations::combiAddSelective(silhouetteLines.size(), allowedCombi);
//		if (print) std::cout << "EEEE combi's: " << combi.size() << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi, silhouetteLines, ray, leaf, edgeRays, printAll, 4, silhouetteTris)) {
//			if (print) std::cout << "EEEE" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}
//
void RaySpaceTree::getEECombis(std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& combi2Edges) {
	if (silhouetteLines.size() >= 2) {
		std::vector<std::vector<int>> combi = Combinations::combi2(silhouetteLines.size());
		for (std::vector<int>& c : combi) {
			Edge e1 = silhouetteEdges[c[0]];
			Edge e2 = silhouetteEdges[c[1]];

			// edge 1 to 2
			bool side1;
			bool check12 = false;
			int silhEdgeCheck = model->checkSilhouetteEdge2(model->verticesIndexed[*e1.vertices.begin()], e2, false, maindir, side1);			
			if (silhEdgeCheck <= 0) {
				bool side2;
				silhEdgeCheck = model->checkSilhouetteEdge2(model->verticesIndexed[*e1.vertices.rbegin()], e2, false, maindir, side2);
				if (silhEdgeCheck == 1);
				else if (silhEdgeCheck == -1) continue;
				else if (side1 == side2) continue;
			}

			// edge 2 to 1
			silhEdgeCheck = model->checkSilhouetteEdge2(model->verticesIndexed[*e2.vertices.begin()], e1, false, maindir, side1);
			if (silhEdgeCheck <= 0) {
				bool side2;
				silhEdgeCheck = model->checkSilhouetteEdge2(model->verticesIndexed[*e2.vertices.rbegin()], e1, false, maindir, side2);
				if (silhEdgeCheck == 1);
				else if (silhEdgeCheck == -1) continue;
				else if (side1 == side2) continue;
			}

			combi2Edges.push_back(c);

		}
	}
}
//
//
//bool RaySpaceTree::checkSEEE(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& splitLines,
//							std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays, 
//							std::vector<std::vector<int>>& allowedCombi) {
//	if (silhouetteLines.size() >= 2 && allowedCombi.size() > 0) {
//		std::vector<std::vector<int>> combi = Combinations::combiAdd1(splitLines.size(), allowedCombi);
//		if (print) std::cout << "SEEE combi's: " << combi.size() << std::endl;
//		if (findExtremalStabbingForPrim(prim, combi, splitLines, ray, leaf, edgeRays, printAll, 3, silhouetteTris, silhouetteLines)) {
//			if (print) std::cout << "SEEE" << std::endl;
//			return true;
//		}
//	}
//	return false;
//}


bool RaySpaceTree::checkEdgeSplittingDuplicates(const int prim,	std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays, std::vector<bool>& sideLines) {
	for (int i = 0; i < sideLines.size(); i++) {
		for (int j = 0; j < edgeRays.size(); j++) {
			if (splitLines[i].equal(edgeRays[j], 1E-5)) { // not very precise?!
				glm::vec3 vx = model->vertices[3 * prim + j].pos;
				glm::vec3 normal_invert = -model->normalPerTri[prim];
				Ray checkRay = Ray(vx + normal_invert, vx);

				bool check1 = splitLines[i].side(checkRay);

				if (splitLines[i].side(checkRay) != sideLines[i]) {
					return false;
				}
			}
		}
	}
	return true;
}


bool RaySpaceTree::checkSilhouetteCombis(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& silhouetteLines,
										std::vector<Ray>& splitLines, std::vector<Edge>& silhouetteEdges, std::vector<int>& silhouetteTris, std::vector<Ray>& triEdgeRays,
										std::set<int>& silhVertices, std::vector<Ray>& silhVertexRays, std::vector<std::vector<int>>& combiV) {


	for (int i = 0; i < silhouetteEdges.size(); i++) {
		silhouetteLines.push_back(Ray(model->verticesIndexed[*silhouetteEdges[i].vertices.begin()], model->verticesIndexed[*silhouetteEdges[i].vertices.rbegin()]));
		for (int t : silhouetteEdges[i].triangles) silhouetteTris.push_back(t);
		if (silhouetteEdges[i].triangles.size() == 1) silhouetteTris.push_back(-1);
		for (int v : silhouetteEdges[i].vertices) {
			if (silhVertices.find(v) == silhVertices.end()) {
				silhVertexRays.push_back(Ray(model->verticesIndexed[v] + model->normalPerTri[silhouetteEdges[i].triangles[0]], model->verticesIndexed[v]));
				combiV.push_back(std::vector<int>{ i, (int) silhVertices.size(), v });
				silhVertices.insert(v);
			}
		}
	}

	std::vector<Ray> e1;
	std::vector<std::vector<int>> e2;
	std::vector<Edge> e3;
	std::vector<int> e4;

	if (checkVeVt(prim, ray, leaf, print, printAll, silhouetteEdges, silhouetteTris, triEdgeRays)) return true;
	
	std::vector<std::vector<int>> combi1S = Combinations::combi1(splitLines.size());

	if (checkCombi(prim, ray, leaf, print, printAll, "SV(E)T", combi1S.size() * combiV.size() * 3, 1, 2, 0, 1,
					splitLines, combi1S, silhouetteLines, e2, e3, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;
	if (checkVeVe(prim, ray, leaf, print, printAll, silhouetteEdges, silhouetteLines, silhouetteTris, triEdgeRays)) return true;
	
	std::vector<std::vector<int>> combi1E = Combinations::combi1(silhouetteLines.size());

	if (checkCombi(prim, ray, leaf, print, printAll, "SV(E)E", combi1S.size() * combi1E.size() * combiV.size(), 1, 2, 1, 0,
					splitLines, combi1S, silhouetteLines, combi1E, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "SEV(T)", combi1S.size() * combi1E.size() * 3, 1, 0, 1, 2,
					splitLines, combi1S, silhouetteLines, combi1E, silhouetteEdges, e1, e2, e4, triEdgeRays)) return true;

	std::vector<std::vector<int>> combi2S = Combinations::combi1(splitLines.size());
	if (checkCombi(prim, ray, leaf, print, printAll, "SSV(E)", combi2S.size() * combiV.size(), 1, 0, 1, 2,
					splitLines, combi2S, silhouetteLines, e2, e3, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "V(E)ET", combi1E.size() * combiV.size() * 3, 0, 2, 1, 1,
					e1, e2, silhouetteLines, combi1E, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;

	std::vector<std::vector<int>> combi2E;
	getEECombis(silhouetteEdges, silhouetteLines, combi2E);
	if (checkCombi(prim, ray, leaf, print, printAll, "EEV(T)", combi2E.size() * 3, 0, 0, 2, 2,
					e1, e2, silhouetteLines, combi1E, silhouetteEdges, e1, e2, e4, triEdgeRays)) return true;

	std::vector<std::vector<int>> combi3E;
	getEEECombis(silhouetteEdges, silhouetteLines, combi2E, combi3E);
	if (checkCombi(prim, ray, leaf, print, printAll, "EEET", combi3E.size() * 3, 0, 0, 3, 1,
					e1, e2, silhouetteLines, combi1E, e3, e1, e2, e4, triEdgeRays)) return true;

	std::vector<std::vector<int>> combi4E = Combinations::combiAddSelective(silhouetteLines.size(), combi3E);
	if (checkCombi(prim, ray, leaf, print, printAll, "EEEE", combi3E.size() * 3, 0, 0, 4, 0,
					e1, e2, silhouetteLines, combi4E, e3, e1, e2, e4, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "V(E)EE", combi2E.size() * combiV.size(), 0, 2, 2, 0,
					e1, e2, silhouetteLines, combi1E, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "SEEE", combi1S.size() * combi3E.size(), 1, 0, 3, 0,
					splitLines, combi1S, silhouetteLines, combi3E, e3, e1, e2, e4, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "SEET", splitLines.size() * combi2E.size() * 3, 1, 0, 2, 1,
					splitLines, combi1S, silhouetteLines, combi1S, silhouetteEdges, e1, e2, e4, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "SSET", combi2S.size() * combi1E.size() * 3, 2, 0, 1, 1,
					splitLines, combi2S, silhouetteLines, combi1E, silhouetteEdges, e1, e2, e4, triEdgeRays)) return true;

	if (checkCombi(prim, ray, leaf, print, printAll, "SSEE", combi2S.size() * combi2E.size(), 2, 0, 2, 0,
					splitLines, combi2S, silhouetteLines, combi2E, e3, e1, e2, e4, triEdgeRays)) return true;

	std::vector<std::vector<int>> combi3S = Combinations::combi3(splitLines.size());
	if (checkCombi(prim, ray, leaf, print, printAll, "SSSE", combi2S.size() * combi2E.size(), 3, 0, 1, 0,
					splitLines, combi2S, silhouetteLines, combi2E, e3, e1, e2, e4, triEdgeRays)) return true;

	return false;
}

bool RaySpaceTree::checkPrim(const int prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, std::vector<std::vector<int>>& combi4,
	std::vector<Ray> splitLines, std::vector<bool>& sideLines, Ray& ray, Node* leaf, bool print, int edgeSelection) {

	glm::vec3 v1 = model->vertices[3 * prim].pos;
	glm::vec3 v2 = model->vertices[3 * prim + 1].pos;
	glm::vec3 v3 = model->vertices[3 * prim + 2].pos;
	std::vector<Ray> edgeRays = { Ray(v3, v2), Ray(v1, v3), Ray(v2, v1) };

	bool printAll = false;
	int splitLinesSize = splitLines.size();

	// Check if triangle edges are same as splitting lines and if yes, if it lies on the correct side
	if (!checkEdgeSplittingDuplicates(prim, splitLines, edgeRays, sideLines)) return false;

	// Check basic combis of extremal stabbing lines
	std::vector<Ray> e1;
	std::vector<std::vector<int>> e2;
	std::vector<Edge> e3;
	std::vector<int> e4;

	//if (checkSSVt(prim, ray, leaf, print, printAll, combi2, splitLines, edgeRays)) return true;
	//if (checkSSST(prim, ray, leaf, print, printAll, combi3, splitLines, edgeRays)) return true;
	if (checkCombi(prim, ray, leaf, print, printAll, "SSV(T)", combi2.size() * 3, 0, 2, 0, 2, splitLines, combi2, e1, e2, e3, e1, e2, e4, edgeRays)) return true;
	if (checkCombi(prim, ray, leaf, print, printAll, "SSST", combi3.size() * 3, 0, 3, 0, 1, splitLines, combi3, e1, e2, e3, e1, e2, e4, edgeRays)) return true;


	// Find silhouette edges for primitive
	std::vector<Edge> silhouetteEdges;
	std::vector<Edge> silhouetteEdgesFirst;
	std::vector<Edge> silhouetteEdgesSecond;

	if (edgeSelection == 0) {
		model->findSilhouetteEdgesForTri(prim, alldir, maindir, silhouetteEdges);
		//exact would be to test if triangles belong in leaf
		for (int i = 0; i < silhouetteEdges.size(); i++) {
			bool found = false;
			for (int pr : silhouetteEdges[i].triangles) {
				if (leaf->primitiveSet.find(pr) != leaf->primitiveSet.end()) {
					silhouetteEdgesFirst.push_back(silhouetteEdges[i]);
					found = true;
					break;
				}
			}
			if (!found) silhouetteEdgesSecond.push_back(silhouetteEdges[i]);
		}
	}
	else if (edgeSelection == 1)
		model->findSilhouetteEdgesForTri(prim, alldir, maindir, silhouetteEdgesFirst, leaf->parent->primitiveSet);

	// Check all combis involving silhouette edges of some sort
	int silhouettesize = 0;
	std::vector<int> silhouetteTris;
	std::vector<Ray> silhouetteLines;
	std::set<int> silhEdgeVertices;
	std::vector<Ray> silhVertexRays;
	std::vector<std::vector<int>> edgeVertexCombis;
	if (checkSilhouetteCombis(prim, ray, leaf, print, printAll, silhouetteLines, splitLines, silhouetteEdgesFirst, 
								silhouetteTris, edgeRays, silhEdgeVertices, silhVertexRays, edgeVertexCombis)) return true;

	// Check basic (but large) combi of extremal stabbing lines
	//if (checkSSSS(prim, ray, leaf, print, printAll, combi4, splitLines, edgeRays)) return true;
	if (checkCombi(prim, ray, leaf, print, printAll, "SSSS", combi4.size(), 0, 4, 0, 0, splitLines, combi4, e1, e2, e3, e1, e2, e4, edgeRays)) return true;

	// Check combis involving second tier silhouette edges
	if (silhouetteEdgesSecond.size() == 0)
		if (checkSilhouetteCombis(prim, ray, leaf, print, printAll, silhouetteLines, splitLines, silhouetteEdgesSecond, 
									silhouetteTris, edgeRays, silhEdgeVertices, silhVertexRays, edgeVertexCombis)) return true;

	if (print) std::cout << "Not Found" << std::endl;
	return false;
}

// Edgeselection 0 --> check from samples
// Edgeselection 1 --> check from parent node
bool RaySpaceTree::checkLeaf(Node* node, std::vector<Ray>& rays, bool getrays, int edgeSelection, bool print) {
	////// TEST TEST /////
	wronglines = std::vector<Ray>();

	std::vector<Ray> splitLines;
	std::vector<bool> sideLines;

	getSplittingLinesInLeafWithSide(node, splitLines, sideLines);
	if (!alldir) {
		std::vector<Ray> boxSides = model->boundingCube.getCubeSideLines(maindir);
		for (auto r : boxSides) splitLines.push_back(r);
	}
	std::vector<Ray> filteredsplitLines;
	std::vector<bool> filterdsideLines;
	filterSplittingLines(node, splitLines, sideLines, filteredsplitLines, filterdsideLines);

	int size = filteredsplitLines.size();// +boxSides.size();
	std::vector<std::vector<int>> combi2 = Combinations::combi2(size);
	std::vector<std::vector<int>> combi3 = Combinations::combi3(size);
	std::vector<std::vector<int>> combi4 = Combinations::combi4(size);

	bool foundAll = true;
	int notfound = 0;

	for (int i : node->primitiveSet) {
		Ray ray;
		if (checkPrim(i, combi2, combi3, combi4, filteredsplitLines, filterdsideLines, ray, node, print, edgeSelection) && getrays) {
			rays.push_back(ray);
		}
		else {
			foundAll = false;
			notfound++;
			std::cout << " DID NOT FIND PRIM " << i << " IN LEAF NR " << node->index << std::endl;

			model->vertices[3 * i].selected = -1.f;
			model->vertices[3 * i + 1].selected = -1.f;
			model->vertices[3 * i + 2].selected = -1.f;
		}
	}
	if (notfound > 0) std::cout << " DID NOT FIND " << notfound << " OF " << node->primitiveSet.size() << " PRIMS IN LEAF NR " << node->index << std::endl;
	model->changeSelected();
	return foundAll;
}


void RaySpaceTree::checkLeaves() {
	int leafcount = 0;
	for (Node* n : nodes) {
		if (!n->leaf) continue;
		std::cout << "leaf index: " << n->index - noLeaves << " from total: " << noLeaves << std::endl;
		std::vector<Ray> rays;
		checkLeaf(n, rays, false, 0);
	}
}

std::vector<Ray> RaySpaceTree::getExtremalStabbingInLeaf(Node* n, bool print) {
	std::vector<Ray> rays;
	checkLeaf(n, rays, true, 0, print);
	for (auto& r : rays) r.get3DfromPlucker();
	return rays;
}