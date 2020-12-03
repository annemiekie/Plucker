#pragma once

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "ray.h"
#include "node.h"
#include "cube.h"
#include <list>
#include "vertex.h"
#include <vector>


class RaySpaceTree {
	Cube* cube;
	Node* rootNode;

public:
	int depth = 0;
	glm::vec3 mainDir = glm::vec3(0);

	std::list<Node> nodes = std::list<Node>();
	std::list<int> visPrims = std::list<int>();

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

	void construct(int level, std::vector<Vertex>& vertices, int option);
	void construct(int level, Node* node, std::vector<Vertex>& vertices, int option);

	std::vector<int> countDuplicates(int size, std::vector<int>& nodenr);
	void countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr);

	void getSplittingLines(Node* node, std::vector<Ray>& lines);
	std::vector<Ray> getSplittingLines();
	std::vector<glm::vec3> getSplittingLinesInCube();

	void getViewingLinesInLeaf(int leafnum, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters);
	void getViewingLinesInLeaf(std::vector<glm::vec3> &rays, std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters, Node* node, int num, int noLeaves, int start);



};