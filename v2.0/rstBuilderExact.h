#pragma once
#include "combinations.h"
#include "rstBuilder.h"
#include "cache.h"
#include "ESLfinder.h"
#include "rst.h"

static bool false_bool = false;

class RSTBuilderExact : public RSTBuilder<RSTBuilderExact> {

public:
	bool cacheEE = false;
	bool cacheEEE = false;
	bool cacheCombi = false;
	Cache<std::vector<Ray>> combiCache;

	static void build(Options::BuildOptions& options, RaySpaceTree* rst, VisComponents& visComp) {//, ESLfinder& eslfinder) {
		
		//rst->model->findPotentialSilhouettes(rst->alldir, rst->maindir);
		//ESLfinder finder(options.cacheCombi);	
		//void RaySpaceTree::fillExact() {
		for (Node* node : rst->nodes) { 
			std::cout << "Computing node nr: " << node->index << " of " << rst->nodes.size() << std::endl;
			if (node->leaf) buildLeaf(rst, node, false, true);
			std::cout << std::endl;
		//	if (cacheCombi)  std::cout << "Combi Cache, Size: " << combiCache.cache.size() << " Hits: " << combiCache.hitcount << " Hash Hit Size: " << combiCache.hashes.size() << std::endl;
		//	if (model->cacheEEE) std::cout << "Edge Edge Edge Cache, Size: " << model->edgeEdgeEdgeCombis.size() << " Hits: " << model->cachehiteee << std::endl;
		//	if (model->cacheEE) std::cout << "Edge Edge Cache, Size: " << model->edgeEdgeCombis[0].size() << " Hits: " << model->cachehitee << std::endl;
		}

	};

	static void fill(RaySpaceTree* rst) {
		for (Node* node : rst->nodes) {
			if (node->leaf) {
				std::cout << "Leaf: " << node->index << std::endl;
				buildLeaf(rst, node, true, true);
				std::cout << std::endl << std::endl;
			}
		}
	};

	static bool find(RaySpaceTree* rst, Line4& ray, Primitive* prim, Node* leaf, bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>()) {
		return check1Prim(rst, prim, ray, leaf, allOptions, allESLs, true);
	}

	static std::vector<Ray> splitLineCombis(RaySpaceTree* rst, Node* leaf);

	static void getSplitLinesInLeaf(RaySpaceTree* rst, Node* leaf, std::vector<Ray>& splitLinesf, std::vector<bool>& sideLinesf,
										bool storecombi = false, std::vector<Ray>& validcombis = std::vector<Ray>());

	static void buildLeaf(RaySpaceTree* rst, Node* leaf, bool fill = true, bool print = false);

	static std::vector<Ray> getExtremalStabbingInLeaf(RaySpaceTree* rst, Node* n, std::vector<int>& notfoundprim = std::vector<int>(), 
														bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(), bool print = false);

	static void checkLeaves(RaySpaceTree* rst);

	static bool checkLeaf(RaySpaceTree* rst, Node* node, std::vector<Ray>& rays, bool getRays, int edgeSelection, 
							std::vector<int>& notfoundprim = std::vector<int>(), bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
							bool print = false);

	static bool check1Prim(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
							bool print = false, int edgeSelection = 0, bool getedges = false,
							std::vector<glm::vec3>& edges = std::vector<glm::vec3>(), std::vector<Ray>& eslEdges = std::vector<Ray>());

	static bool checkPrim(RaySpaceTree* rst, Primitive* prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, 
							std::vector<std::vector<int>>& combi4, std::vector<Ray> splitLines, std::vector<bool>& sideLines, Line4& ray, Node* leaf, 
							bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(), bool print = false, int edgeSelection = 0,
							bool getedges = false, std::vector<glm::vec3>& edges = std::vector<glm::vec3>(), std::vector<Ray>& eslEdges = std::vector<Ray>());

	static void filterSplittingLines(RaySpaceTree* rst, Node* leaf, std::vector<Ray>& splitlines, std::vector<bool>& sides,
										std::vector<Ray>& filteredLines, std::vector<bool>& filteredSides, bool storecombi = false,
										std::vector<Ray>& validcombis = std::vector<Ray>());

	static bool checkExtremalStabbingLine(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, 
											bool print = false, bool& inBox = false_bool,
											std::vector<Ray>& lines = std::vector<Ray>(), bool vischeck = true,
											std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(),
											std::vector<Edge*>& silhEdges = std::vector<Edge*>());

	static bool checkRayAndReverse(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print = false,
									bool& inBox = false_bool, std::vector<Ray>& lines = std::vector<Ray>(), bool vischeck = true,
									std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(),
									std::vector<Edge*>& silhEdges = std::vector<Edge*>());


	static bool checkRaysThroughLines(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print = false,
										std::vector<Ray>& lines = std::vector<Ray>(), std::vector<uint64_t>& indices = std::vector<uint64_t>(),
										bool vischeck = true, std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(),
										std::vector<Edge*>& silhEdges = std::vector<Edge*>());

	static bool checkRayInPrim(Line4& ray, std::vector<Ray>& lines4, Primitive* prim, bool& inPlane, bool print);

	static bool checkRayInPlane(Line4& ray, Primitive* prim, bool& inPlane, bool print);

	static bool checkRayInBox(RaySpaceTree* rst, const Line4& ray, bool print = false);

	static bool checkRayInLeaf(RaySpaceTree* rst, Node* node, const Line4& ray, std::vector<Ray>& lines, bool print = false);

	static bool checkNoSplitVertexProblem(RaySpaceTree* rst, Node* node, Primitive *prim, const Line4& ray, std::vector<Ray>& lines, bool print = false);

	static bool checkPrimVisibleForRay(RaySpaceTree* rst, Line4& ray, Primitive* prim, std::vector<Vertex*>& silhVertices, 
										std::vector<Edge*>& silhEdges, std::vector<Ray>& lines, bool inplane, bool print);

	static bool checkEdgeSplittingDuplicates(Primitive* prim, std::vector<Ray>& splitLines, std::vector<bool>& sideLines);

	static bool checkEdgeInLeafCombis(RaySpaceTree* rst, Edge* e, Node* leaf, std::vector<Ray>& splitLines, 
										 std::vector<std::vector<int>>& combi, Line4& ray);

	static bool checkVeVt(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Vertex*>& silhVertices,
							bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>());

	static bool checkVeVe(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Vertex*>& silhVertices,
							bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>());

	static bool checkSilhouetteCombis(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll,
										std::vector<Ray>& splitLines, std::vector<Edge*>& silhouetteEdgesToAdd, 
										std::vector<Edge*>& silhouetteEdges, std::set<Vertex*, Vertex::cmp_ptr>& silhVertices, 
										bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>());

	static bool checkCombi(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, std::string combi_text, int combiNr,
							int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges, 
							std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis, bool vischeck = true, 
							bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
							std::vector<std::vector<int>>& silhLineCombis = std::vector<std::vector<int>>(), 
							std::vector<Edge*>& silhouetteEdges = std::vector<Edge*>(),
							std::vector<Vertex*> silhVertices = std::vector<Vertex*>());

	static void getEECombis(RaySpaceTree* rst, std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges);
	static void getEEECombis(RaySpaceTree* rst, std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges,
								std::vector<std::vector<int>>& combi3Edges);

	static bool checkSilhouettesForRay(Line4& ray, std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges,
										 std::map<float, float>& intersectionDepthOffset);

};
