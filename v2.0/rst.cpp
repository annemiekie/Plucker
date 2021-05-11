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
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * (model->vertices.size() / 3));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int redge = int(r2 * 3.f);
		float r3 = (float)rand() / static_cast <float> (RAND_MAX);
		float roffset = 0.1 * r3;
		splitter = Ray(model->vertices[3 * r1 + redge].pos + roffset, model->vertices[3 * r1 + (redge + 1) % 3].pos);
	}
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

void RaySpaceTree::constructSmartRandom(glm::ivec2 res, std::vector<Camera*>& cams,
										std::vector<std::pair<int, int>>& samples, std::set<int>& tris) {
	int t = (int)time(NULL);
	//srand(t);
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
	//nodeSamples next = { tris, samples, rootNode, 0 };
	//toProcess.push(next);

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
	//if (abs(splitter.sideVal(ray)) < 1E-8) std::cout << "OH NOOOOOO " << abs(splitter.sideVal(ray)) <<" ";
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

void RaySpaceTree::putPrimitive(Ray& ray, int primId, bool putRay) {
	putPrimitive(ray, primId, rootNode, putRay);
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, Node* node, bool putRay) {
	if (node->leaf) {
		if (putRay) node->insert(primId, ray);
		else node->insert(primId);
		//visPrims.insert(primId);
		return;
	}
	if (node->splitter.side(ray)) putPrimitive(ray, primId, node->leftNode, putRay);
	else putPrimitive(ray, primId, node->rightNode, putRay);
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

//void RaySpaceTree::descendWithLines(Ray& ray, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, GeoObject* object)
//{
//	return descendWithLines(ray, rootNode, splitters, rays, colors, object);
//}
//
//void RaySpaceTree::descendWithLines(Ray& ray, Node* node, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, GeoObject* object)
//{
//	if (node->leaf) return getViewingLinesInLeafInGeo(node, rays, colors, object);
//	addLineInGeo(node->splitter, object, splitters);
//	if (node->splitter.side(ray)) return descendWithLines(ray, node->leftNode, splitters, rays, colors, object);
//	else return descendWithLines(ray, node->rightNode, splitters, rays, colors, object);
//}

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

//void RaySpaceTree::addSplittersForLeaf(std::vector<glm::vec3>& splitters, Node* node, GeoObject* object) {
//	addSplitter(splitters, node, object);
//	if (node != rootNode) return addSplittersForLeaf(splitters, node->parent, object);
//	return;
//}
//
//void RaySpaceTree::addSplitter(std::vector<glm::vec3>& splitters, Node* node, GeoObject* object) {
//	glm::vec3 s, t;
//	object->intersect(node->splitter, s, t);
//
//	splitters.push_back(s);
//	splitters.push_back(t);
//}

// add ignore list
bool RaySpaceTree::checkLineInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool print) {

	for (auto er : edgeRays) {
		bool equal = false;
		for (auto igray : lines4) {
			if (igray.equal(er, 1E-8)) { //-8
				equal = true;
				break;
			}
		}
		if (!equal && er.sideVal(line) < -1E-8) { //-8
			if (print) std::cout << std::endl << "On the wrong side: " << er.sideVal(line) << std::endl;
			return false;
		}
	}
	glm::dvec3 cross = glm::cross(glm::normalize(edgeRays[0].direction), glm::normalize(edgeRays[1].direction));
	line.get3DfromPlucker();
	double v = fabsf(glm::dot(line.direction, cross));
	if (v < 1E-10) { //-10
		if (print) std::cout << std::endl << "In plane: " << v << std::endl;
		for (int i = 0; i < 3; i++) {
			glm::dvec3 vert = model->vertices[3 * prim + i].pos;
			glm::dvec3 dir = vert - line.origin;
			if (dir.x * line.direction.x < 0) dir = -dir;
			double dist = glm::distance(glm::normalize(dir), glm::normalize(line.direction));
			if (print) std::cout << "Through vertex? " << dist << std::endl;
			if (dist < 1E-10) return true; //-10
		}
		return false;
	}
	return true;
	//return !(v < 1E-10);
}


bool RaySpaceTree::rayInLeaf(Node* node, Ray& ray, std::vector<Ray>& ignore, bool print) {
	if (node == rootNode) return true;
	Node* parent = node->parent;
	bool ign = false;
	for (auto igray : ignore) {
		if (igray.equal(parent->splitter, 1E-8)) ign = true; //-8
	}
	double sidev = ray.sideVal(parent->splitter);
	if (abs(sidev) < 1E-8) ign = true; //-15
	
	if (print) std::cout << std::endl << "RayInLeaf for node: " << node->index << " value: " << sidev << " ignore: " << ign << std::endl;
	//bool side = ray.side(parent->splitter);
	if (parent->leftNode == node) {
		if (sidev < 0 || ign) return rayInLeaf(parent, ray, ignore, print);
		else return false;
	}
	else {
		if (sidev > 0 || ign) return rayInLeaf(parent, ray, ignore, print);
		else return false;
	}
}

//bool RaySpaceTree::checkLineInBox(Ray& ray) {
//	//float planeVal = glm::dot(glm::abs(mainDir), cube->bounds[b])
//	return model->boundingCube.intersectSide(mainDir, ray);
//}

bool RaySpaceTree::checkPrim(const int prim, std::vector<std::vector<int>>& splitCombi, std::vector<Ray>& splitLines, Ray &ray, Node* leaf) {
	bool print = false;
	bool printAll = false;
	bool details = false;
	//if (leaf->index == 41) print == true;
	//if (prim == 42 && print) printAll = true;

	glm::vec3 v1 = model->vertices[3 * prim].pos;
	glm::vec3 v2 = model->vertices[3 * prim + 1].pos;
	glm::vec3 v3 = model->vertices[3 * prim + 2].pos;
	std::vector<Ray> edgeRays = { Ray(v2, v1), Ray(v3, v2), Ray(v1, v3) };
	std::vector<std::vector<int>> edges = { {2,0}, {0,1}, {1,2} };
	int edgeCombi = 3;
	if (splitCombi[0].size() == 4) edgeCombi = 1;

	bool found = false;
	for (int e = 0; e < edgeCombi; e++) {
		if (found) break;
		for (int c = 0; c < splitCombi.size(); c++) {
			if (found) break;
			// this check can be done earlier for speedup
			bool dupli = false;
			for (int i = 0; i < splitCombi[c].size(); i++) {
				if (dupli) break;
				for (int j = 0; j<4-splitCombi[c].size(); j++) {
					if (splitLines[splitCombi[c][i]].equal(edgeRays[edges[e][j]], 1E-10)) {
						dupli = true;
						break;
					}
				}
			}
			if (dupli) continue;

			std::vector<Ray> lines;
			std::vector<Ray> ignore;
			if (splitCombi[0].size() == 2) {
				lines = { edgeRays[edges[e][0]], edgeRays[edges[e][1]], splitLines[splitCombi[c][0]], splitLines[splitCombi[c][1]] };
				ignore = { splitLines[splitCombi[c][0]], splitLines[splitCombi[c][1]] };
			}
			else if (splitCombi[0].size() == 3) {
				lines = { edgeRays[e], splitLines[splitCombi[c][0]], splitLines[splitCombi[c][1]], splitLines[splitCombi[c][2]] };
				ignore = { splitLines[splitCombi[c][0]], splitLines[splitCombi[c][1]], splitLines[splitCombi[c][2]] };
			}
			else if (splitCombi[0].size() == 4) {
				lines = { splitLines[splitCombi[c][0]], splitLines[splitCombi[c][1]], splitLines[splitCombi[c][2]], splitLines[splitCombi[c][3]] };
				ignore = lines;
			}

			std::vector<Ray> intersectLines = LineThroughFour::find(lines, model);

			for (int i = 0; i < intersectLines.size() * 2; i++) {
				if (found) break;
				Ray line = intersectLines[i/2];
				if (i % 2 == 1) line.inverseDir(); // check both directions of extremal stabbing line
				if (printAll) {
					std::cout << "Leaf:" << leaf->index << " Prim:" << prim << " c:" << c << " edges: ";
					edgeRays[edges[e][0]].print();
					std::cout << " and ";
					edgeRays[edges[e][1]].print(); 
					std::cout << std::endl << " -- Testing ray: ";
					line.print();
					std::cout << " test in prim: " << checkLineInPrim(edgeRays, line, lines, prim, printAll && (c == 5 || c == 7));
					std::cout << " test in leaf: " << rayInLeaf(leaf, line, ignore, printAll && (c==5 || c==7));
					std::cout << std::endl;
				}
				//if (!checkLineInBox(line)) continue;
				if (!checkLineInPrim(edgeRays, line, lines, prim, false)) continue;
				if (!rayInLeaf(leaf, line, ignore, details)) continue;
				ray = line;
				found = true;
				if (print) {
					std::cout << "Leaf:" << leaf->index << " Prim:" << prim << " c:" << c << " edges: ";
					edgeRays[edges[e][0]].print();
					std::cout << " and ";
					edgeRays[edges[e][1]].print();
					std::cout << std::endl << " found ray ";
					ray.print();
					std::cout << std::endl;
				}
			}
		}
	}
	return found;
}

bool RaySpaceTree::checkLeaf(Node* node, std::vector<Ray>& rays, bool getrays) {
	std::vector<Ray> splitLines;
	getSplittingLinesInLeaf(node, splitLines);
	//std::vector<Ray> boxSides = model->boundingBox.getCubeSideLines(mainDir);
	//for (auto r : boxSides) splitLines.push_back(r);
	int size = splitLines.size();// +boxSides.size();
	std::vector<std::vector<int>> combi2 = Combinations::combi2(size);
	std::vector<std::vector<int>> combi3 = Combinations::combi3(size);
	std::vector<std::vector<int>> combi4 = Combinations::combi4(size);
	

	bool foundAll = true;
	int notfound = 0;

	//for (int i = 0; i < model->vertices.size(); i++) {
	//	model->vertices[i].selected = 0.f;
	//}

	for (int i : node->primitiveSet) {
		Ray ray;
		bool found = checkPrim(i, combi2, splitLines, ray, node);
		if (!found) found = checkPrim(i, combi3, splitLines, ray, node); // only for occlusion!!!
	//	if (!found) found = checkPrim(i, combi4, splitLines, ray, node); // only for occlusion!!!
		if (!found) {
			foundAll = false;
			notfound++;
			std::cout << " DID NOT FIND PRIM " << i << " IN LEAF NR " << node->index << std::endl;

			model->vertices[3 * i].selected = -1.f;
			model->vertices[3 * i + 1].selected = -1.f;
			model->vertices[3 * i + 2].selected = -1.f;
		}
		if (getrays && found)  {
			rays.push_back(ray);
		}
	}
	if (notfound > 0) std::cout << " DID NOT FIND " << notfound << " OF " << node->primitiveSet.size() << " PRIMS IN LEAF NR " << node->index << std::endl;
	model->changeSelected();
	return foundAll;
}

void RaySpaceTree::checkLeaves() {
	int leafcount = 0;
	for (Node *n : nodes) {
		if (!n->leaf) continue;
		std::cout << "leaf index: " << n->index << " from total: " << noLeaves << std::endl;
		std::vector<Ray> rays;
		checkLeaf(n, rays, false);
	}
}

std::vector<Ray> RaySpaceTree::getExtremalStabbingInLeaf(Node* n) {
	std::vector<Ray> rays;
	checkLeaf(n, rays, true);
	for (auto &r : rays) r.get3DfromPlucker();
	return rays;
}
