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

	static bool find(RaySpaceTree* rst, Line4& ray, Primitive* prim, Node* leaf,
		bool allESLs = false, std::vector<Line4>& esls = std::vector<Line4>()) {
		return check1Prim(rst, prim, ray, leaf, allESLs, esls, true);
	};

	static bool check1Prim(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf,
		bool allESLs, std::vector<Line4>& allESLsForPrimNode, bool print) {//, int edgeSelection,
					//bool getedges, std::vector<glm::vec3>& edges, std::vector<Ray>& esledges) {

		std::vector<Split> splitLines = getSplitLinesInLeaf(rst, leaf);
		if (splitLines.size() == 0) return;
		CombiConfigurations combiS(splitLines.size());
		ESLFinder eslFinder(rst, prim, leaf, false, print, allESLs, splitLines, combiS, false);

		if (eslFinder.find()) {
			if (allESLs) allESLsForPrimNode = eslFinder.allESLsForPrimNode;
			ray = eslFinder.ray;
			return true;
		}
		else if (print) std::cout << "Not Found" << std::endl;
		return false;
		//if (checkPrim(rst, prim, combi2, combi3, combi4, splitLines, sideLines, ray, leaf,
		//	allOptions, allESLs, print, edgeSelection, getedges, edges, esledges) || (allOptions && allESLs.size())) return true;


	};

	//static bool checkLeaf(RaySpaceTree* rst, Node* leaf, std::vector<Ray>& rays, bool getrays, bool checkAllPrimsForSilh,
	//						std::vector<int>& notfoundprim = std::vector<int>(), 
	//						bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
	//						bool print = false) {
	static bool checkLeaf(RaySpaceTree* rst, Node* leaf, std::vector<Ray>& rays, bool getrays, bool checkAllPrimsForSilh,
							std::vector<int>& notfoundprim = std::vector<int>(),
							bool allESLs = false, std::vector<Ray>& allESLsForPrimNode = std::vector<Ray>(), bool print = false) {

		std::vector<Split> splitLines = getSplitLinesInLeaf(rst, leaf);
		int size = splitLines.size();
		if (size == 0) {
			if (leaf->primitiveSet.size() == 0) return true;
			else return false;
		}
		CombiConfigurations combiS(size);

		for (int id : leaf->primitiveSet) {
			ESLFinder eslFinder(rst, rst->model->triangles[id], leaf, false, print, allESLs, splitLines, combiS, false);
			if (eslFinder.find() && getrays) rays.push_back(eslFinder.ray);
			else {
				notfoundprim.push_back(id);
				std::cout << " DID NOT FIND PRIM " << id << " IN LEAF NR " << leaf->index << std::endl;
			}
		}
		if (notfoundprim.size() > 0) std::cout << " DID NOT FIND " << notfoundprim.size() << " OF "
												<< leaf->primitiveSet.size() << " PRIMS IN LEAF NR "
												<< leaf->index << std::endl;
		return notfoundprim.size() == 0;
	}

	// check where to put the "only regard primitives in parent nodes" ///
	static void buildLeaf(RaySpaceTree* rst, Node* leaf, bool fill, bool print) {
		std::vector<Split> splitLines = getSplitLinesInLeaf(rst, leaf);
		if (splitLines.size() == 0) return;
		CombiConfigurations combiS(splitLines.size());

		for (Primitive* prim : rst->model->triangles) {
			if (print) std::cout << prim->id << ", ";
			if (fill && leaf->primitiveSet.find(prim->id) != leaf->primitiveSet.end()) continue;
			ESLFinder eslFinder(rst, prim, leaf, false, print, false, splitLines, combiS, false);
			if (eslFinder.find()) leaf->insert(prim->id);
			//if (checkPrim(rst, prim, combi2, combi3, combi4, splitLines, sideLines, r, leaf)) leaf->insert(prim->id);
		}
	}

	//static std::vector<Ray> splitLineCombis(RaySpaceTree* rst, Node* leaf);

	//static void filterSplittingLines(RaySpaceTree* rst, Node* leaf, std::vector<Ray>& splitlines, std::vector<bool>& sides,
	//	std::vector<Ray>& filteredLines, std::vector<bool>& filteredSides, bool storecombi = false,
	//	std::vector<Ray>& validcombis = std::vector<Ray>());
	static std::vector<Split> filterSplittingLines(RaySpaceTree* rst, Node* leaf, std::vector<Split>& splitlines,
													bool storecombi = false, std::vector<Ray>& validCombis = std::vector<Ray>()) {

		std::vector<int> foundExtremalStabbing;
		std::vector<std::vector<int>>& splitCombi4 = Combinations::combi4(splitlines.size());
		ESLChecker eslcheck(leaf, NULL, rst);

		for (int i = 0; i < splitCombi4.size(); i++) {
			std::vector<Ray> lines4;
			for (int c : splitCombi4[i]) lines4.push_back(splitlines[c].line);
			std::vector<Line4> intersectLines = Lines4Finder::find(lines4);
			for (Line4& r : intersectLines) {
				ESLCandidate esl({ r, lines4 });
				if (eslcheck.isExtremalStabbing(esl, false, false)) {
					if (storecombi) validCombis.push_back(r);
					foundExtremalStabbing.push_back(i);
					break;
				}
			}
		}
		std::vector<int> splitLineTimes(splitlines.size());
		std::set<int> filtered;
		for (int esl : foundExtremalStabbing) {
			for (int i : splitCombi4[esl]) {
				filtered.insert(i);
				splitLineTimes[i]++;
			}
		}
		std::vector<Split> filteredLines;
		//for (int i : splitLineTimes) if (i == foundExtremalStabbing.size()) return;
		for (int i : filtered) filteredLines.push_back(splitlines[i]);
			//if (i < splitlines.size() - 4) filteredSides.push_back(sides[i]);
		return filteredLines;
	}

	static std::vector<Split> getSplitLinesInLeaf(RaySpaceTree* rst, Node* leaf,
													bool storeCombi = false, std::vector<Ray>& validCombis = std::vector<Ray>()) {

		std::vector<Split> splitLines;

		rst->getSplittingLinesInNodeWithSide(leaf, splitLines);
		if (!rst->alldir) {
			std::vector<Ray> boxSides = rst->model->boundingCube.getCubeSideSquare(rst->maindir).quadLines;
			for (int i = 0; i < 4; i++) {
				Ray boxside = boxSides[i];
				boxside.index = splitLines.size() + i;
				splitLines.push_back(Split({ boxside }));
			}
		}

		return filterSplittingLines(rst, leaf, splitLines, storeCombi, validCombis);
	}



	static std::vector<Ray> splitLineCombis(RaySpaceTree* rst, Node* leaf) {
		std::vector<Ray> validCombis;
		getSplitLinesInLeaf(rst, leaf, true, validCombis);
		for (Ray& r : validCombis) r.get3DfromPlucker();
		return validCombis;
	}



	static void RSTBuilderExact::checkLeaves(RaySpaceTree* rst) {
		int leafcount = 0;
		for (Node* n : rst->nodes) {
			if (!n->leaf) continue;
			std::cout << "leaf index: " << n->index - rst->noLeaves << " from total: " << rst->noLeaves << std::endl;
			std::vector<Ray> rays;
			checkLeaf(rst, n, rays, false, true);
		}
	}

	static std::vector<Ray> getExtremalStabbingInLeaf(RaySpaceTree* rst, Node* n, std::vector<int>& notfoundprim,
														bool allOptions, std::vector<Ray>& allESLs, bool print) {
		std::vector<Ray> rays;
		checkLeaf(rst, n, rays, true, 0, notfoundprim, allOptions, allESLs, print);
		for (auto& r : rays) r.get3DfromPlucker();
		for (auto& r : allESLs) r.get3DfromPlucker();
		return rays;
	}

	//static void buildLeaf(RaySpaceTree* rst, Node* leaf, bool fill = true, bool print = false);

	//static std::vector<Ray> getExtremalStabbingInLeaf(RaySpaceTree* rst, Node* n, std::vector<int>& notfoundprim = std::vector<int>(), 
	//													bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(), bool print = false);

	//static void checkLeaves(RaySpaceTree* rst);

	//static bool checkLeaf(RaySpaceTree* rst, Node* node, std::vector<Ray>& rays, bool getRays, int edgeSelection, 
	//						std::vector<int>& notfoundprim = std::vector<int>(), bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
	//						bool print = false);

	//static bool check1Prim(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
	//						bool print = false, int edgeSelection = 0, bool getedges = false,
	//						std::vector<glm::vec3>& edges = std::vector<glm::vec3>(), std::vector<Ray>& eslEdges = std::vector<Ray>());

	//static bool checkPrim(RaySpaceTree* rst, Primitive* prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, 
	//						std::vector<std::vector<int>>& combi4, std::vector<Ray> splitLines, std::vector<bool>& sideLines, Line4& ray, Node* leaf, 
	//						bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(), bool print = false, int edgeSelection = 0,
	//						bool getedges = false, std::vector<glm::vec3>& edges = std::vector<glm::vec3>(), std::vector<Ray>& eslEdges = std::vector<Ray>());



	//static bool checkExtremalStabbingLine(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, 
	//										bool print = false, bool& inBox = false_bool,
	//										std::vector<Ray>& lines = std::vector<Ray>(), bool vischeck = true,
	//										std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(),
	//										std::vector<Edge*>& silhEdges = std::vector<Edge*>());

	//static bool checkRayAndReverse(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print = false,
	//								bool& inBox = false_bool, std::vector<Ray>& lines = std::vector<Ray>(), bool vischeck = true,
	//								std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(),
	//								std::vector<Edge*>& silhEdges = std::vector<Edge*>());


	//static bool checkRaysThroughLines(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print = false,
	//									std::vector<Ray>& lines = std::vector<Ray>(), std::vector<uint64_t>& indices = std::vector<uint64_t>(),
	//									bool vischeck = true, std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(),
	//									std::vector<Edge*>& silhEdges = std::vector<Edge*>());

	//static bool checkRayInPrim(Line4& ray, std::vector<Ray>& lines4, Primitive* prim, bool& inPlane, bool print);

	//static bool checkRayInPlane(Line4& ray, Primitive* prim, bool& inPlane, bool print);

	//static bool checkRayInBox(RaySpaceTree* rst, const Line4& ray, bool print = false);

	//static bool checkRayInLeaf(RaySpaceTree* rst, Node* node, const Line4& ray, std::vector<Ray>& lines, bool print = false);

	//static bool checkNoSplitVertexProblem(RaySpaceTree* rst, Node* node, Primitive *prim, const Line4& ray, std::vector<Ray>& lines, bool print = false);

	//static bool checkPrimVisibleForRay(RaySpaceTree* rst, Line4& ray, Primitive* prim, std::vector<Vertex*>& silhVertices, 
	//									std::vector<Edge*>& silhEdges, std::vector<Ray>& lines, bool inplane, bool print);

	//static bool checkEdgeSplittingDuplicates(Primitive* prim, std::vector<Ray>& splitLines, std::vector<bool>& sideLines);

	//static bool checkEdgeInLeafCombis(RaySpaceTree* rst, Edge* e, Node* leaf, std::vector<Ray>& splitLines, 
	//									 std::vector<std::vector<int>>& combi, Line4& ray);

	//static bool checkVeVt(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Vertex*>& silhVertices,
	//						bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>());

	//static bool checkVeVe(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Vertex*>& silhVertices,
	//						bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>());

	//static bool checkSilhouetteCombis(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll,
	//									std::vector<Ray>& splitLines, std::vector<Edge*>& silhouetteEdgesToAdd, 
	//									std::vector<Edge*>& silhouetteEdges, std::set<Vertex*, Vertex::cmp_ptr>& silhVertices, 
	//									bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>());

	//static bool checkCombi(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, std::string combi_text, int combiNr,
	//						int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges, 
	//						std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis, bool vischeck = true, 
	//						bool allOptions = false, std::vector<Ray>& allESLs = std::vector<Ray>(),
	//						std::vector<std::vector<int>>& silhLineCombis = std::vector<std::vector<int>>(), 
	//						std::vector<Edge*>& silhouetteEdges = std::vector<Edge*>(),
	//						std::vector<Vertex*> silhVertices = std::vector<Vertex*>());

	//static void getEECombis(RaySpaceTree* rst, std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges);
	//static void getEEECombis(RaySpaceTree* rst, std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges,
	//							std::vector<std::vector<int>>& combi3Edges);

	//static bool checkSilhouettesForRay(Line4& ray, std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges,
	//									std::set<float>& intersectionDepths);

};
