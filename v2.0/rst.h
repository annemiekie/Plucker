#pragma once

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "ray.h"
#include "node.h"
#include "cube.h"
#include "vertex.h"
#include <vector>
#include "model.h"
#include "nodeSamples.h"
#include "orthocamera.h"
#include "model.h"
#include "combinations.h"
#include "textureRenderer.h"
#include <embree3/rtcore.h>
#include "cache.h"

class RaySpaceTree {

public:
	//Cache<std::vector<Line4>> combiCache;
	bool cacheCombi = false;
	Node* rootNode;
	int depth = 0;
	int noLeaves = 0;
	glm::vec3 maindir;
	bool alldir = false;
	std::vector<Node*> nodes = std::vector<Node*>();
	std::set<int> visPrims = std::set<int>();
	std::vector<Ray> wronglines = std::vector<Ray>();
	Model* model;
	std::vector<Ray> splitters;

	RaySpaceTree() {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};

	RaySpaceTree(Model* model, int depth, bool alldir = false, glm::vec3 maindir = glm::vec3())
		: model(model), depth(depth), alldir(alldir), maindir(maindir) {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};
	~RaySpaceTree() {};

	void putPrimitive(Ray& ray, int primId, bool putRay, bool putPrim = true);
	void putPrimitive(Ray& ray, int primId, Node* node, bool putRay, bool putPrim = true);

	Node* descend(Ray& ray);
	Node* descend(Ray& ray, Node* node);

	void printLeafNodes();

	bool intoLeftNode(Ray& splitter, std::vector<Camera*>& cams, int raynr, glm::ivec2 res);
	Ray getRay(std::vector<Camera*>& cams, int raynr, glm::ivec2 res);

	std::vector<int> countDuplicates(int size, std::vector<int>& nodenr);
	void countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr);

	void getSplittingLines(Node* node, std::vector<Ray>& lines);
	std::vector<Ray> getAllSplittingLines();
	std::vector<glm::vec3> getSplittingLinesInGeo(GeoObject* object);

	std::vector<Ray> getViewingLinesInLeaf(int leafnum);
	std::vector<Ray> getViewingLinesInLeaf(Node* n);
	std::vector<Ray> getviewingLinesInLeafInTri(Node* n, int prim);
	Node* getLeafFromNum(int leafNum);

	void getSplittingLinesInLeaf(int leafnum, std::vector<Ray>& lines);
	void getSplittingLinesInNode(Node* n, std::vector<Ray>& lines);
	void getSplittingLinesInNodeWithSide(Node* n, std::vector<Ray>& lines, std::vector<bool>& sides);

	int numOfLeaf(int ind);
	int getNumberOfTriInleaf(int leafnum);

	void printTree();
	void printTree(Node* node, int level);

};