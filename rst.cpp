#include "rst.h"

void RaySpaceTree::construct(int lvl, std::vector<Vertex>& vertices, int option) {
	int t = (int)time(NULL);
	srand(t);
	construct(lvl, rootNode, vertices, option);
}

void RaySpaceTree::construct(int lvl, Node *node, std::vector<Vertex>& vertices, int option) {
	//rootNode.splitter = Ray::getRandom();
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


void RaySpaceTree::putPrimitive(Ray& ray, int primId) {
	putPrimitive(ray, primId, rootNode);
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, Node* node) {
	if (node->leaf) {
		node->insert(primId, ray);
		visPrims.push_back(primId);
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
		glm::vec3 s, t, x;
		cube->intersect(node->splitter, s, t, mainDir, x);
		splitters.push_back(s);
		splitters.push_back(t);
		if (num < start + noLeaves / 2) getViewingLinesInLeaf(rays, colors, splitters, node->leftNode, num, noLeaves / 2, start);
		else getViewingLinesInLeaf(rays, colors, splitters, node->rightNode, num, noLeaves / 2, start + noLeaves/2);
	}
	else {
		for (std::pair<int, Ray> p : node->primAndRayVector) {
			glm::vec3 start, end, planeEnter;
			cube->intersect(p.second, start, end, mainDir, planeEnter);
			rays.push_back(start);
			rays.push_back(end);
			colors.push_back(planeEnter);
		}
	}
}

