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
		rst->filledExact = true;
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
		rst->filledExact = true;
		for (Node* node : rst->nodes) {
			if (node->leaf) {
				std::cout << "Leaf: " << node->index << std::endl;
				buildLeaf(rst, node, true, true);
				std::cout << std::endl << std::endl;
			}
		}
	};

	static bool check1Prim(RaySpaceTree* rst, Primitive* prim, Node* leaf,
							bool allESLs = false, std::vector<Ray>& esls = std::vector<Ray>(),
							bool print = false) {//, int edgeSelection,

		std::vector<Split> splitLines = getSplitLinesInLeaf(rst, leaf);
		if (splitLines.size() == 0) return false;
		CombiConfigurations combiS(splitLines.size());
		ESLFinder eslFinder(rst, prim, leaf, false, print, allESLs, splitLines, combiS, !rst->filledExact);

		if (eslFinder.find()) {
			for (Ray& r : eslFinder.esls) esls.push_back(r);
			return true;
		}
		else if (print) std::cout << "Not Found" << std::endl;
		return false;
	};

	static bool checkLeaf(RaySpaceTree* rst, Node* leaf, std::vector<Ray>& rays, bool getrays, 
							std::vector<int>& notfoundprim = std::vector<int>(),
							bool storeAllEsls = false, bool print = false) {

		std::vector<Split> splitLines = getSplitLinesInLeaf(rst, leaf);
		int size = splitLines.size();
		if (size == 0) {
			if (leaf->primitiveSet.size() == 0) return true;
			else return false;
		}
		CombiConfigurations combiS(size);

		for (int id : leaf->primitiveSet) {
			ESLFinder eslFinder(rst, rst->model->triangles[id], leaf, false, print, storeAllEsls, splitLines, combiS, !rst->filledExact);
			if (!eslFinder.find()) {
				notfoundprim.push_back(id);
				std::cout << " DID NOT FIND PRIM " << id << " IN LEAF NR " << leaf->index << std::endl;
			}			
			else if (getrays) for (Ray& r : eslFinder.esls) rays.push_back(r);
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
		}
	}

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
			checkLeaf(rst, n, rays, false);
		}
	}

	static std::vector<Ray> getExtremalStabbingInLeaf(RaySpaceTree* rst, Node* n, std::vector<int>& notfoundprim,
														bool storeAllESLs, bool print = false) {
		std::vector<Ray> rays;
		checkLeaf(rst, n, rays, true, notfoundprim, storeAllESLs, print);
		return rays;
	}
};