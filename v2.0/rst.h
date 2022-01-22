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
#include "model.h"
#include <queue>
#include "nodeSamples.h"
#include "orthocamera.h"
#include "model.h"
#include "findLineThroughFour.h"
#include "combinations.h"
#include "textureRenderer.h"
#include <embree3/rtcore.h>
#include "cache.h"

static bool false_bool = false;

class RaySpaceTree {
	std::queue<nodeSamples> toProcess = std::queue<nodeSamples>();

public:
	Cache<std::vector<Line4>> combiCache;
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

	void getStatistics();

	//void descendWithLines(Ray& ray, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, GeoObject* object);
	//void descendWithLines(Ray& ray, Node* node, std::vector<glm::vec3>& splitters, std::vector<glm::vec3>& rays, std::vector<glm::vec3>& colors, GeoObject* object);

	void construct(int option, std::vector<Ray>& rays);
	void construct(int level, Node* node, int option, std::vector<Ray>& rays, std::vector<Ray> splitters, int splitnum);

	void fillExact();
	void fillExact(Node* node);

	void constructAdaptive(glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples, 
							std::set<int>& tris, bool print);
	void constructAdaptive(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples, 
							std::set<int>& tris, bool print, std::vector<Ray> splitters);
	void constructAdaptive2(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples, 
							std::set<int>& tris, int totSampSize);

	void constructSmartRandom(glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples, std::set<int>& tris);
	void constructSmartRandom(int level, Node* node, glm::ivec2 res, std::vector<Camera*>& cams, std::vector<std::pair<int, int>>& samples, 
								std::set<int>& tris, std::vector<Ray> splitters);

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
	std::vector<Ray>  getViewingLinesInLeaf(Node* n);
	Node* getLeafFromNum(int leafNum);

	void getSplittingLinesInLeaf(int leafnum, std::vector<Ray>& lines);
	void getSplittingLinesInLeaf(Node* n, std::vector<Ray>& lines);
	void getSplittingLinesInLeafWithSide(Node* n, std::vector<Ray>& lines, std::vector<bool>& sides);

	int numOfLeaf(int ind);
	int getNumberOfTriInleaf(int leafnum);
	//uint64_t makeCombiKey(std::vector<Ray>& lines, int nrOfSplitLines);
	//uint64_t makeCombiKey(std::vector<Ray>& lines, int nrOfSplitLines, int nrOfSilhVertices);

	std::vector<Ray> getExtremalStabbingInLeaf(Node* n, std::vector<int>& notfoundprim = std::vector<int>(), bool print = false);
	void checkLeaves();
	bool checkLeaf(Node* node, std::vector<Ray>& rays, bool getRays, int edgeSelection, std::vector<int>& notfoundprim = std::vector<int>(), bool print = false);
	bool check1Prim(const int prim, Ray& ray, Node* leaf, bool print, int edgeSelection, bool getedges = false, std::vector<glm::vec3>& edges = std::vector<glm::vec3>(),
					std::vector<Ray>& eslEdges = std::vector<Ray>());
	bool checkPrim(const int prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, std::vector<std::vector<int>>& combi4,
					std::vector<Ray> splitLines, std::vector<bool>& sideLines, Ray& ray, Node* leaf, bool print, int edgeSelection,
					bool getedges = false, std::vector<glm::vec3>& edges = std::vector<glm::vec3>(), std::vector<Ray>& eslEdges = std::vector<Ray>());

	void filterSplittingLines(Node* leaf, std::vector<Ray>& splitlines, std::vector<bool>& sides,
		std::vector<Ray>& filteredLines, std::vector<bool>& filteredSides);

	bool checkExtremalStabbingLine(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
									bool print = false, std::vector<int>& visibleTriIgnore = std::vector<int>(), bool& inBox = false_bool,
									std::vector<Ray>& lines = std::vector<Ray>(), int rayIgnoresize = 0, bool vischeck = true);
	bool checkRayAndReverse(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
							bool print = false, std::vector<int>& visibleTriIgnore = std::vector<int>(), bool& inBox = false_bool,
							std::vector<Ray>& lines = std::vector<Ray>(), int rayIgnoresize = 0, bool vischeck = true);
	bool checkRaysThroughLines(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
								bool print = false, std::vector<int>& visibleTriIgnore = std::vector<int>(),
								std::vector<Ray>& lines = std::vector<Ray>(), std::vector<uint64_t>& indices = std::vector<uint64_t>(),
								int rayIgnoresize = 0, bool vischeck = true);
	bool checkRayInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool& inPlane, bool print);
	bool checkRayInBox(const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print);
	bool checkRayInLeaf(Node* node, const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print);

	bool checkPrimVisibleForRay(Ray& ray, const int prim, std::vector<int>& ignore, std::vector<Ray>& lines, bool inplane, bool print);

	bool checkEdgeSplittingDuplicates(const int prim, std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays, std::vector<bool>& sideLines);

	bool checkEdgeInLeafCombis(Edge& e, Node* leaf, std::vector<Ray>& splitLines, std::string combi_text,
								int nrOfsplitLines, std::vector<std::vector<int>>& combi, Ray& ray);

	bool checkVeVt(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges,
					std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays);

	bool checkVeVe(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges,
					std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays);

	bool checkSilhouetteCombis(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& silhouetteLines,
								std::vector<Ray>& splitLines, std::vector<Edge>& silhouetteEdgesToAdd, std::vector<Edge>& silhouetteEdges, std::vector<int>& silhouetteTris,
								std::vector<Ray>& edgeRays, std::set<int>& silhEdgeVertices, std::vector<Ray>& edgeVertexRays,
								std::vector<std::vector<int>>& edgeVertexCombis);

	bool checkCombi(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::string combi_text, int combiNr,
					int nrOfsilhEdges, int nrOfsplitLines, int nrOfVertices, int nrOfTriEdges,
					std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis,
					std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& silhLineCombis, std::vector<Edge>& silhouetteEdges,
					std::vector<Ray>& silhVertexLines, std::vector<std::vector<int>>& silhVertexCombis,
					std::vector<int>& silhouetteTris, std::vector<Ray>& triEdgeRays, bool vischeck = true);

	void getEECombis(std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& combi2Edges);
	void getEEECombis(std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& combi2Edges, 
						std::vector<std::vector<int>>& combi3Edges);

	void printTree();
	void printTree(Node* node, int level);

};