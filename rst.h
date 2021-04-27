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
#include <queue>
#include "nodeSamples.h"
#include "lineThroughFour.h"
#include "combinations.h"

class RaySpaceTree {


	std::queue<nodeSamples> toProcess = std::queue<nodeSamples>();

public:
	Cube* cube;

	Node* rootNode;
	Model* model;
	int depth = 0;
	int noLeaves = 0;
	glm::vec3 mainDir = glm::vec3(0);

	std::vector<Node*> nodes = std::vector<Node*>();
	std::set<int> visPrims = std::set<int>();

	RaySpaceTree() {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};

	RaySpaceTree(Cube* cube) : cube(cube) {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};

	RaySpaceTree(int depth, Cube* cube, glm::vec3 mainDir) : depth(depth), mainDir(mainDir), cube(cube) {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};
	~RaySpaceTree() {};

	void putPrimitive(Ray& ray, int primId, bool putRay);
	void putPrimitive(Ray& ray, int primId, Node* node, bool putRay);

	Node* descend(Ray& ray);
	Node* descend(Ray& ray, Node* node);

	void descendWithLines(Ray& ray, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors);
	void descendWithLines(Ray& ray, Node* node, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors);

	void construct(int option, std::vector<Ray>& rays);
	void construct(int level, Node* node, int option, std::vector<Ray>& rays, std::vector<Ray> splitters);
	
	void constructAdaptive(glm::ivec2 res, std::vector<Orthocamera>& cams,
							std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print);
	void constructAdaptive(int level, Node* node, glm::ivec2 res, std::vector<Orthocamera>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray> splitters);
	void constructAdaptive2(int level, Node* node, glm::ivec2 res, std::vector<Orthocamera>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize);

	void constructSmartRandom(glm::ivec2 res, std::vector<Orthocamera>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris);
	void constructSmartRandom(int level, Node* node, glm::ivec2 res, std::vector<Orthocamera>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, std::vector<Ray> splitters);

	bool getBestSplitter(Ray& splitter, glm::ivec2 res, std::vector<std::pair<int, int>>& samplesL, std::set<int>& trisL,
		std::vector<std::pair<int, int>>& samplesR, std::set<int>& trisR, std::vector<Orthocamera>& cams, int tries,
						std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray>& splitters);
	Ray getBestSplitter2(glm::ivec2 res, std::vector<Orthocamera>& cams, int tries,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize);

	bool intoLeftNode(Ray& splitter, std::vector<Orthocamera>& cams, int raynr, glm::ivec2 res);
	Ray getRay(std::vector<Orthocamera>& cams, int raynr, glm::ivec2 res);

	std::vector<int> countDuplicates(int size, std::vector<int>& nodenr);
	void countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr);

	void getSplittingLines(Node* node, std::vector<Ray>& lines);
	std::vector<Ray> getSplittingLines();
	std::vector<glm::vec3> getSplittingLinesInCube(Cube* cube);

	void getViewingLinesInLeaf(int leafnum, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters, Cube* cube);
	//void getViewingLinesInLeaf(std::vector<glm::vec3> &rays, std::vector<glm::vec3>& colors, std::vector<glm::vec3>& splitters, Node* node, int num, int start);
	Node* getLeafFromNum(int leafNum);
	void getSplittingLinesInLeaf(int leafnum, std::vector<Ray>& lines);
	void getSplittingLinesInLeaf(Node* n, std::vector<Ray>& lines);

	void addSplitter(std::vector<glm::vec3>& splitters, Node* node, Cube* cube);
	void addSplittersForLeaf(std::vector<glm::vec3>& splitters, Node* node, Cube* cube);

	int numOfLeaf(int ind);

	void addViewingLines(std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, Node* node);

	void checkLeaves();
	bool checkLeaf(Node *node, std::vector<Ray>& rays, bool getRays);
	bool checkPrim(int i, std::vector<std::vector<int>>& splitCombi, std::vector<Ray>& splitLines, Ray &ray, Node *node);
	bool checkLineInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool print);
	bool checkLineInBox(Ray& ray);
	bool rayInLeaf(Node* node, Ray& ray, std::vector<Ray>& ignore, bool print);

	void printTree();
	void printTree(Node* node, int level);

	void postProcessNeighbours();

};