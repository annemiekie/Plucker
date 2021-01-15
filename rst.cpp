#include "rst.h"
#include <chrono>

void RaySpaceTree::construct(std::vector<Vertex>& vertices, int option) {
	int t = (int)time(NULL);
	srand(t);
	construct(0, rootNode, vertices, option);
}

void RaySpaceTree::construct(int lvl, Node *node, std::vector<Vertex>& vertices, int option) {
	if (lvl >= depth) {
		node->leaf = true;
		return;
	}
	if (option == 0) { // edges
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * (vertices.size() / 3));
		r1 *= 3;
		node->splitter = Ray(vertices[r1].pos, vertices[r1 + 1].pos);
	}
	else if (option == 1) { // ortogonal from vertex
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * vertices.size());
		glm::vec3 pos = vertices[r1].pos;
		if (rand() % 2) node->splitter = Ray(glm::vec3(pos.x, 0, pos.z), glm::vec3(pos.x, 1, pos.z));
		else node->splitter = Ray(glm::vec3(pos.x, pos.y, 0), glm::vec3(pos.x, pos.y, 1));
	}
	else if (option == 2) { // ortogonal random
		float rX = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		float rZ = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		float rY = 2.f * float(rand()) / static_cast <float> (RAND_MAX);

		if (rand() % 2) node->splitter = Ray(glm::vec3(rX, 0, rZ), glm::vec3(rX, 1, rZ));
		else node->splitter = Ray(glm::vec3(rX, rY, 0), glm::vec3(rX, rY, 1));
	}


	nodes.push_back(Node(node->index*2 + 1));
	node->leftNode = &(nodes.back());
	construct(lvl + 1, node->leftNode, vertices, option);

	nodes.push_back(Node(node->index*2 + 2));
	node->rightNode = &(nodes.back());
	construct(lvl + 1, node->rightNode, vertices, option);
}

void RaySpaceTree::constructAdaptive(std::vector<Vertex>& vertices, glm::ivec2 res, std::vector<Orthocamera>& cams, 
	std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize, int zerosize)
{
	constructAdaptive(0, rootNode, vertices, res, cams, samples, tris, totSampSize, zerosize);
}

void RaySpaceTree::constructAdaptive(int level, Node* node, std::vector<Vertex>& vertices, glm::ivec2 res, std::vector<Orthocamera>& cams,
									std::vector<std::pair<int,int>>& samples, std::set<int>& tris, int totSampSize, int zerosize) {
	auto start_time = std::chrono::high_resolution_clock::now();
	std::cout << "level: " << level << std::endl;
	if (level < depth) {
		std::cout << "Node with " << samples.size() << " samples" << std::endl;
		std::cout << "Finding best splitter..." << std::endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		Ray splitter;
		if (getBestSplitter(splitter, vertices, res, cams, 10, samples, tris, totSampSize, zerosize)) {

			auto end_time = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
			std::cout << "Found splitter in " << diff << " ms" << std::endl;

			int zerosizeL = 0;
			int zerosizeR = 0;
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
					if (tri >= 0) {
						trisL.insert(tri);
					}
					else {
						zerosizeL++;
					}
				}
				else {
					samplesR.push_back(sample);
					if (tri >= 0) {
						trisR.insert(tri);
					}
					else { 
						zerosizeR++;
					}
				}
			}
			end_time = std::chrono::high_resolution_clock::now();
			diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
			std::cout << "Samples put left and right in " << diff << " ms" << std::endl;


			nodes.push_back(Node(node->index * 2 + 1));
			node->leftNode = &(nodes.back());
			constructAdaptive(level + 1, node->leftNode, vertices, res, cams, samplesL, trisL, totSampSize, zerosizeL);

			nodes.push_back(Node(node->index * 2 + 2));
			node->rightNode = &(nodes.back());
			constructAdaptive(level + 1, node->rightNode, vertices, res, cams, samplesR, trisR, totSampSize, zerosizeR);

			return;
		}
	}
	std::cout << "Leaf with " << samples.size() << " samples";
	node->leaf = true;
	node->primitiveSet = tris;
	for (auto sample : samples) {
		if (sample.second >= 0) node->insert(sample.second, getRay(cams, sample.first, res));
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	std::cout << "Filled leaf in " << diff << " ms" << std::endl;
}

Ray RaySpaceTree::getRay(std::vector<Orthocamera>& cams, int raynr, glm::ivec2 res) {
	int camNr= raynr / (res.x*res.y);

	int rayInCam = raynr % (res.x * res.y);
	int y = rayInCam / res.x;
	int x = rayInCam % res.x;
	const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
	return cams[camNr].pixRayDirection(pixpos);
}

bool RaySpaceTree::intoLeftNode(Ray &splitter, std::vector<Orthocamera>& cams, int raynr, glm::ivec2 res) {
	Ray ray = getRay(cams, raynr, res);
	return splitter.side(ray);
}

bool RaySpaceTree::getBestSplitter(Ray& splitter, std::vector<Vertex>& vertices, glm::ivec2 res, std::vector<Orthocamera>& cams, int tries,
									std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize, int zeroSize) {
	int t = (int)time(NULL);
	srand(t);
	// NOPE: Number of samples that do not hit everything are counted twice, and penalized more that way.
	float fraction = float(samples.size()) / float(totSampSize);
	float bestSAH = tris.size() * fraction;
	std::cout << "Parent SAH = " << bestSAH << std::endl;

	for (int i = 0; i < tries; i++) {
		std::cout << "Checking line no " << i << "... ";
		// do not test splitting lines that already have been tested???
		// ..needs to be implemented

		// choose the first egde of a random triangle in this node
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int rtri = int(r * tris.size());
		int redge = int(r * 3);
		std::set<int>::iterator it = tris.begin();
		std::advance(it, rtri);
		int tri = *it;
		Ray testRay = Ray(vertices[3*tri + redge].pos, vertices[3* tri + (redge + 1) % 3].pos);

		//int r2 = int(r * (vertices.size() / 3));
		//Ray test = Ray(vertices[3 * r2].pos, vertices[3 * r2 + 1].pos);

		std::set<int> trisLeft = std::set<int>();
		std::set<int> trisRight = std::set<int>();
		int pixLeft = 0;
		int pixRight = 0;
		//int zeroLeft = 0;
		//int zeroRight = 0;

		for (auto sample : samples) {
			bool left = intoLeftNode(testRay, cams, sample.first, res);
			int tri = sample.second;
			if (left) {
				pixLeft++;
				if (tri >= 0) trisLeft.insert(tri);
				//else zeroLeft++;
			}
			else {
				pixRight++;
				if (tri >= 0) trisRight.insert(tri);
				//else zeroRight++;
			}
		}
		float fracLeft = 1.f * (pixLeft) / float(totSampSize);
		float fracRight = 1.f * (pixRight) / float(totSampSize);

		float sahleft = trisLeft.size() * fracLeft; // sah
		float sahright = trisRight.size() * fracRight;
		std::cout << " for a SAH of " << sahleft+sahright;

		if (sahleft + sahright < bestSAH) {
			std::cout << " which is better!" << std::endl;
			bestSAH = sahleft + sahright;
			splitter = testRay;
		}
		else {
			std::cout << std::endl;
		}
	}
	std::cout << "Old value: " << tris.size() * fraction << ", New value: " << bestSAH << std::endl;
	if (tris.size() * fraction == bestSAH) return false;
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

void RaySpaceTree::putPrimitive(Ray& ray, int primId) {
	putPrimitive(ray, primId, rootNode);
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, Node* node) {
	if (node->leaf) {
		node->insert(primId, ray);
		//visPrims.insert(primId);
		return;
	}
	if (node->splitter.side(ray)) putPrimitive(ray, primId, node->leftNode);
	else putPrimitive(ray, primId, node->rightNode);
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

void RaySpaceTree::descendWithLines(Ray& ray, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors)
{
	return descendWithLines(ray, rootNode, splitters, rays, colors);
}

void RaySpaceTree::descendWithLines(Ray& ray, Node* node, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors)
{
	if (node->leaf) return addViewingLines(rays, colors, node);
	addSplitter(splitters, node);
	if (node->splitter.side(ray)) return descendWithLines(ray, node->leftNode, splitters, rays, colors);
	else return descendWithLines(ray, node->rightNode, splitters, rays, colors);
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

std::vector<Ray> RaySpaceTree::getSplittingLines() {
	std::vector<Ray> lines = std::vector<Ray>();
	getSplittingLines(rootNode, lines);
	return lines;
}

void RaySpaceTree::getSplittingLines(Node* node, std::vector<Ray>& lines) {
	if (!node->leaf) {
		lines.push_back(node->splitter);
		getSplittingLines(node->leftNode, lines);
		getSplittingLines(node->rightNode, lines);
	}
}

std::vector<glm::vec3> RaySpaceTree::getSplittingLinesInCube() {
	std::vector<Ray> lines = getSplittingLines();
	std::vector<glm::vec3> linePos;
	glm::vec3 start, end, x;
	for (Ray line : lines) {
		cube->intersect(line, start, end, mainDir, x);
		linePos.push_back(start);
		linePos.push_back(end);
	}
	return linePos;
}

void RaySpaceTree::getViewingLinesInLeaf(int leafNum, std::vector<glm::vec3>& rays, 
										std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters)
{
	int noLeaves = std::pow(2, depth);
	getViewingLinesInLeaf(rays, colors, splitters, rootNode,leafNum, noLeaves, 0);
}

void RaySpaceTree::getViewingLinesInLeaf(std::vector<glm::vec3> &rays, std::vector<glm::vec3>& colors,
										std::vector<glm::vec3>& splitters, Node *node, int num, int noLeaves, int start)
{
	if (!node->leaf) {
		addSplitter(splitters, node);
		if (num < start + noLeaves / 2) getViewingLinesInLeaf(rays, colors, splitters, node->leftNode, num, noLeaves / 2, start);
		else getViewingLinesInLeaf(rays, colors, splitters, node->rightNode, num, noLeaves / 2, start + noLeaves/2);
	}
	else addViewingLines(rays, colors, node);
}

void RaySpaceTree::addSplitter(std::vector<glm::vec3>& splitters, Node* node) {
	glm::vec3 s, t, x;
	cube->intersect(node->splitter, s, t, mainDir, x);
	splitters.push_back(s);
	splitters.push_back(t);
}

void RaySpaceTree::addViewingLines(std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, Node* node) {
	for (Sample p : node->primAndRayVector) {
		glm::vec3 start, end, planeEnter;
		cube->intersect(p.ray, start, end, mainDir, planeEnter);
		rays.push_back(start);
		rays.push_back(end);
		colors.push_back(planeEnter);
	}
}