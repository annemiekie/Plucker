//#pragma once
//#include "findLineThroughFour.h"
//#include "combinations.h"
//#include "rstBuilder.h"
//#include "cache.h"
//#include "node.h"
//static bool false_bool = false;
//
//class RSTBuilderExact : public RSTBuilder<RSTBuilderExact>{
//
//public:
//	bool cacheEE = false;
//	bool cacheEEE = false;
//	bool cacheCombi = false;
//	Cache<std::vector<Ray>> combiCache;
//
//	void build(BuildOptions& options, RaySpaceTree* rst) {};
//
//
//	//void fillExact(RaySpaceTree *rst);
//	//void fillExact(Node* node);
//
//	//std::vector<Ray> getExtremalStabbingInLeaf(Node* n, std::vector<int>& notfoundprim = std::vector<int>(), bool print = false);
//
//	//void checkLeaves();
//
//	//bool checkLeaf(Node* node, std::vector<Ray>& rays, bool getRays, int edgeSelection, std::vector<int>& notfoundprim = std::vector<int>(), bool print = false);
//
//	//bool check1Prim(const int prim, Ray& ray, Node* leaf, bool print, int edgeSelection, bool getedges = false, std::vector<glm::vec3>& edges = std::vector<glm::vec3>(),
//	//	std::vector<Ray>& eslEdges = std::vector<Ray>());
//
//	//bool checkPrim(const int prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, std::vector<std::vector<int>>& combi4,
//	//	std::vector<Ray> splitLines, std::vector<bool>& sideLines, Ray& ray, Node* leaf, bool print, int edgeSelection,
//	//	bool getedges = false, std::vector<glm::vec3>& edges = std::vector<glm::vec3>(), std::vector<Ray>& eslEdges = std::vector<Ray>());
//
//	//void filterSplittingLines(Node* leaf, std::vector<Ray>& splitlines, std::vector<bool>& sides,
//	//	std::vector<Ray>& filteredLines, std::vector<bool>& filteredSides);
//
//	//bool checkExtremalStabbingLine(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
//	//	bool print = false, std::vector<int>& visibleTriIgnore = std::vector<int>(), bool& inBox = false_bool,
//	//	std::vector<Ray>& lines = std::vector<Ray>(), int rayIgnoresize = 0, bool vischeck = true);
//
//	//bool checkRayAndReverse(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
//	//	bool print = false, std::vector<int>& visibleTriIgnore = std::vector<int>(), bool& inBox = false_bool,
//	//	std::vector<Ray>& lines = std::vector<Ray>(), int rayIgnoresize = 0, bool vischeck = true);
//
//	//bool checkRaysThroughLines(const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
//	//	bool print = false, std::vector<int>& visibleTriIgnore = std::vector<int>(),
//	//	std::vector<Ray>& lines = std::vector<Ray>(), std::vector<uint64_t>& indices = std::vector<uint64_t>(),
//	//	int rayIgnoresize = 0, bool vischeck = true);
//
//	//bool checkRayInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool& inPlane, bool print);
//
//	//bool checkRayInBox(const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print);
//
//	//bool checkRayInLeaf(Node* node, const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print);
//
//	//bool checkPrimVisibleForRay(Ray& ray, const int prim, std::vector<int>& ignore, std::vector<Ray>& lines, bool inplane, bool print);
//
//	//bool checkEdgeSplittingDuplicates(const int prim, std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays, std::vector<bool>& sideLines);
//
//	//bool checkEdgeInLeafCombis(Edge& e, Node* leaf, std::vector<Ray>& splitLines, std::string combi_text,
//	//	int nrOfsplitLines, std::vector<std::vector<int>>& combi, Ray& ray);
//
//	//bool checkVeVt(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges,
//	//	std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays);
//
//	//bool checkVeVe(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges,
//	//	std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays);
//
//	//bool checkSilhouetteCombis(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& silhouetteLines,
//	//	std::vector<Ray>& splitLines, std::vector<Edge>& silhouetteEdgesToAdd, std::vector<Edge>& silhouetteEdges, std::vector<int>& silhouetteTris,
//	//	std::vector<Ray>& edgeRays, std::set<int>& silhEdgeVertices, std::vector<Ray>& edgeVertexRays,
//	//	std::vector<std::vector<int>>& edgeVertexCombis);
//
//	//bool checkCombi(const int prim, Ray& ray, Node* leaf, bool print, bool printAll, std::string combi_text, int combiNr,
//	//	int nrOfsilhEdges, int nrOfsplitLines, int nrOfVertices, int nrOfTriEdges,
//	//	std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis,
//	//	std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& silhLineCombis, std::vector<Edge>& silhouetteEdges,
//	//	std::vector<Ray>& silhVertexLines, std::vector<std::vector<int>>& silhVertexCombis,
//	//	std::vector<int>& silhouetteTris, std::vector<Ray>& triEdgeRays, bool vischeck = true);
//
//	//void getEECombis(std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& combi2Edges);
//	//void getEEECombis(std::vector<Edge>& silhouetteEdges, std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& combi2Edges,
//	//	std::vector<std::vector<int>>& combi3Edges);
//
//};
