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

	Node* rootNode;
	Model* model;
	int depth = 0;
	int noLeaves = 0;
	glm::vec3 maindir;

	std::vector<Node*> nodes = std::vector<Node*>();
	std::set<int> visPrims = std::set<int>();

	RaySpaceTree() {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};

	RaySpaceTree(int depth) : depth(depth) {
		rootNode = new Node(0, 0);
		nodes.push_back(rootNode);
	};
	~RaySpaceTree() {};

	void putPrimitive(Ray& ray, int primId, bool putRay);
	void putPrimitive(Ray& ray, int primId, Node* node, bool putRay);

	Node* descend(Ray& ray);
	Node* descend(Ray& ray, Node* node);

	//void descendWithLines(Ray& ray, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, GeoObject* object);
	//void descendWithLines(Ray& ray, Node* node, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, GeoObject* object);

	void construct(int option, std::vector<Ray>& rays);
	void construct(int level, Node* node, int option, std::vector<Ray>& rays, std::vector<Ray> splitters);
	
	void constructAdaptive(glm::ivec2 res, std::vector<Camera*>& cams,
							std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print);
	void constructAdaptive(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray> splitters);
	void constructAdaptive2(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize);

	void constructSmartRandom(glm::ivec2 res, std::vector<Camera*>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris);
	void constructSmartRandom(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, std::vector<Ray> splitters);

	bool getBestSplitter(Ray& splitter, glm::ivec2 res, std::vector<std::pair<int, int>>& samplesL, std::set<int>& trisL,
		std::vector<std::pair<int, int>>& samplesR, std::set<int>& trisR, std::vector<Camera*>& cams, int tries,
						std::vector<std::pair<int, int>>& samples, std::set<int>& tris, bool print, std::vector<Ray>& splitters);
	Ray getBestSplitter2(glm::ivec2 res, std::vector<Camera*>& cams, int tries,
		std::vector<std::pair<int, int>>& samples, std::set<int>& tris, int totSampSize);

	bool intoLeftNode(Ray& splitter, std::vector<Camera*>& cams, int raynr, glm::ivec2 res);
	Ray getRay(std::vector<Camera*>& cams, int raynr, glm::ivec2 res);

	std::vector<int> countDuplicates(int size, std::vector<int>& nodenr);
	void countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr);

	void getSplittingLines(Node* node, std::vector<Ray>& lines);
	std::vector<Ray> getAllSplittingLines();
	std::vector<glm::vec3> getSplittingLinesInGeo(GeoObject* object);

	std::vector<Ray>  getViewingLinesInLeaf(int leafnum);
	std::vector<Ray>  getViewingLinesInLeaf(Node *n);
	Node* getLeafFromNum(int leafNum);
	void getSplittingLinesInLeaf(int leafnum, std::vector<Ray>& lines);
	void getSplittingLinesInLeaf(Node* n, std::vector<Ray>& lines);
	//void addSplitter(std::vector<glm::vec3>& splitters, Node* node, GeoObject* object);
	//void addSplittersForLeaf(std::vector<glm::vec3>& splitters, Node* node, GeoObject* object);

	int numOfLeaf(int ind);
	int getNumberOfTriInleaf(int leafnum);

	//void addViewingLines(std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, Node* node, GeoObject* object);

	std::vector<Ray> getExtremalStabbingInLeaf(Node* n);
	void checkLeaves();
	bool checkLeaf(Node *node, std::vector<Ray>& rays, bool getRays);
	bool checkPrim(int i, std::vector<std::vector<int>>& splitCombi, std::vector<Ray>& splitLines, Ray &ray, Node *node);
	bool checkLineInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool print);
	//bool checkLineInBox(Ray& ray);
	bool rayInLeaf(Node* node, Ray& ray, std::vector<Ray>& ignore, bool print);

	void printTree();
	void printTree(Node* node, int level);

	//void postProcessNeighbours();

};