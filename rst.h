#pragma once

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "ray.h"
#include "node.h"
#include "cube.h"
#include <list>
#include <set>
#include "vertex.h"
#include <vector>
#include "orthocamera.h"
#include "model.h"


class RaySpaceTree {
	Cube* cube;

public:
	Node* rootNode;

	int depth = 0;
	glm::vec3 mainDir = glm::vec3(0);

	std::list<Node> nodes = std::list<Node>();
	std::set<int> visPrims = std::set<int>();

	RaySpaceTree() {
		nodes.push_back(Node(0));
		rootNode = &(nodes.back());		
	};

	RaySpaceTree(Cube* cube) : cube(cube) {
		nodes.push_back(Node(0));
		rootNode = &(nodes.back());
	};

	RaySpaceTree(int depth, Cube* cube, glm::vec3 mainDir) : depth(depth), mainDir(mainDir), cube(cube) {
		nodes.push_back(Node(0));
		rootNode = &(nodes.back());
	};
	~RaySpaceTree() {};

	void putPrimitive(Ray& ray, int primId);
	void putPrimitive(Ray& ray, int primId, Node* node);

	Node* descend(Ray& ray);
	Node* descend(Ray& ray, Node* node);

	void descendWithLines(Ray& ray, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors);
	void descendWithLines(Ray& ray, Node* node, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors);

	void construct(std::vector<Vertex>& vertices, int option);
	void construct(int level, Node* node, std::vector<Vertex>& vertices, int option);

	void constructAdaptive(std::vector<Vertex>& vertices, glm::ivec2 res, std::vector<Orthocamera>& cams, 
							std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize, int zerosize);
	void constructAdaptive(int level, Node* node, std::vector<Vertex>& vertices, glm::ivec2 res, std::vector<Orthocamera>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize, int zerosize);

	bool getBestSplitter(Ray& splitter, std::vector<Vertex>& vertices, glm::ivec2 res, std::vector<Orthocamera>& cams, int tries,
						std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize, int zerosize);
	bool intoLeftNode(Ray& splitter, std::vector<Orthocamera>& cams, int raynr, glm::ivec2 res);
	Ray getRay(std::vector<Orthocamera>& cams, int raynr, glm::ivec2 res);

	std::vector<int> countDuplicates(int size, std::vector<int>& nodenr);
	void countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr);

	void getSplittingLines(Node* node, std::vector<Ray>& lines);
	std::vector<Ray> getSplittingLines();
	std::vector<glm::vec3> getSplittingLinesInCube();

	void getViewingLinesInLeaf(int leafnum, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters);
	void getViewingLinesInLeaf(std::vector<glm::vec3> &rays, std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters, Node* node, int num, int noLeaves, int start);

	void addSplitter(std::vector<glm::vec3>& splitters, Node* node);
	void addViewingLines(std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, Node* node);

	void printTree();
	void printTree(Node* node, int level);

};