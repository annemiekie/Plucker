#pragma once
#include "combinations.h"
#include "rstBuilder.h"
#include "cache.h"
#include "ESLfinder.h"
#include "rst.h"
#include <chrono>
#include "exactBuildTracker.h"


static bool false_bool = false;

class RSTBuilderExact : public RSTBuilder<RSTBuilderExact> {

public:
	//bool cacheEE = false;
	//bool cacheEEE = false;
	//bool cacheCombi = false;
	//Cache<std::vector<Ray>> combiCache;

	static void build(Options::BuildOptions& options, RaySpaceTree* rst, VisComponents& visComp, bool print) {
		fill(options, rst, print);
	};

	static void fill(Options::BuildOptions& options, RaySpaceTree* rst, bool print, bool depthFirst = false) {
		ExactBuildTracker track;
		std::deque<Node*> nodesToProcess;
		nodesToProcess.push_back(rst->rootNode);

		while (!nodesToProcess.empty()) {
			Node* node;
			if (depthFirst) node = nodesToProcess.back();
			else node = nodesToProcess.front();

			if (node->depth <= rst->depth) {
				if (node->depth >= options.exactStartLevel) {
					if (print) std::cout 
						<< "Filling Node: " << node->index << " / " << rst->nodes.size() << ". Progress: ";
					auto start_time = std::chrono::system_clock::now();
					fillNode(rst, node, print, track);
					auto end_time = std::chrono::system_clock::now();
					int time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
					track.totalTime += time / 1000;
					if (print) std::cout << std::endl 
						<< "Time to fill node: " << time << " ms, Total time: " << track.totalTime << " s"
						<< std::endl << std::endl;
				}
				if (!node->leaf) {
					nodesToProcess.push_back(node->leftNode);
					nodesToProcess.push_back(node->rightNode);
				}
			}
			if (depthFirst) nodesToProcess.pop_back();
			else nodesToProcess.pop_front();
		}

		//	//	if (cacheCombi)  std::cout << "Combi Cache, Size: " << combiCache.cache.size() << " Hits: " << combiCache.hitcount << " Hash Hit Size: " << combiCache.hashes.size() << std::endl;
		//	//	if (model->cacheEEE) std::cout << "Edge Edge Edge Cache, Size: " << model->edgeEdgeEdgeCombis.size() << " Hits: " << model->cachehiteee << std::endl;
		//	//	if (model->cacheEE) std::cout << "Edge Edge Cache, Size: " << model->edgeEdgeCombis[0].size() << " Hits: " << model->cachehitee << std::endl;
		rst->filledExact = true;
	}


	//static void fill(RaySpaceTree* rst, bool print = false) {
	//
	//	for (Node* node : rst->nodes) {
	//		if (node->leaf) {
	//			std::cout << "Leaf: " << node->index << std::endl;
	//			fillNode(rst, node, true, print);
	//			std::cout << std::endl << std::endl;
	//		}
	//	}
	//	rst->filledExact = true;
	//};

	static bool check1Prim(RaySpaceTree* rst, Primitive* prim, Node* leaf,
							bool allESLs = false, std::vector<Ray>& esls = std::vector<Ray>(),
							bool print = false) {//, int edgeSelection,

		std::vector<SplitSide> splitLines = getSplitLinesInNode(rst, leaf);
		if (splitLines.size() == 0) return false;
		CombiConfigurations combiS(splitLines.size());
		ESLFinder eslFinder(rst, prim, leaf, false, print, allESLs, splitLines, combiS);

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

		std::vector<SplitSide> splitLines = getSplitLinesInNode(rst, leaf);
		int size = splitLines.size();
		if (size == 0) {
			if (leaf->primitiveSet.size() == 0) return true;
			else return false;
		}
		CombiConfigurations combiS(size);

		for (int id : leaf->primitiveSet) {
			ESLFinder eslFinder(rst, rst->model->triangles[id], leaf, false, print, storeAllEsls, splitLines, combiS);
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
	static void fillNode(RaySpaceTree* rst, Node* node, bool print, ExactBuildTracker& track) {

		std::vector<SplitSide> splitLines = getSplitLinesInNode(rst, node);
		if (splitLines.size() == 0) return;
		CombiConfigurations combiS(splitLines.size());
		std::vector<Primitive*> checkPrims;
		if (node->parent->filledExact) {
			for (int primId : node->parent->primitiveSet) checkPrims.push_back(rst->model->triangles[primId]);
		}
		else checkPrims = rst->model->triangles;

		float nrOfSteps = 10;
		float percentagePerStep = 100.f / nrOfSteps;
		float primPerStep = checkPrims.size() / nrOfSteps;
		int primnum = 0;
		int step = 0;

		for (Primitive* prim : checkPrims) {
			if (print) {
				if (primnum - int(step * primPerStep) == 0) {
					std::cout << int(percentagePerStep * step) << "%...";
					step++;
				}
				primnum++;
			}

			if (node->primitiveSet.find(prim->id) != node->primitiveSet.end()) {
				track.alreadyInNode++;
				continue;
			}
			ESLFinder eslFinder(rst, prim, node, false, false, false, splitLines, combiS);
			if (eslFinder.find()) {
				if (eslFinder.checkHardCombis) track.inNodeSlow++;
				else track.inNodeFast++;
				node->insert(prim->id);
			}
			else {
				if (eslFinder.checkHardCombis) track.notInNodeSlow++;
				else track.notInNodeFast++;
			}
		}
		track.totalPrimitives += checkPrims.size();
		if (print) std::cout << "100%" << std::endl;
		if (print) std::cout << "Total primitives checked: " << track.totalPrimitives
								<< ", Already in node: " << track.alreadyInNode
								<< ", Not in node fast: " << track.notInNodeFast
								<< ", slow: " << track.notInNodeSlow
								<< ", In node fast: " << track.inNodeFast
								<< ", slow: " << track.inNodeSlow;


		node->filledExact = true;
	}

	static std::vector<SplitSide> filterSplittingLines(RaySpaceTree* rst, Node* leaf, std::vector<SplitSide>& splitlines,
													bool storecombi = false, std::vector<Ray>& validCombis = std::vector<Ray>()) {

		std::vector<int> foundExtremalStabbing;
		std::vector<std::vector<int>>& splitCombi4 = Combinations::combi4(splitlines.size());
		ESLChecker eslcheck(leaf, NULL, rst);

		for (int i = 0; i < splitCombi4.size(); i++) {
			std::vector<Ray> lines4;
			for (int c : splitCombi4[i]) lines4.push_back(splitlines[c].ray);
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
		std::vector<SplitSide> filteredLines;
		//for (int i : splitLineTimes) if (i == foundExtremalStabbing.size()) return;
		for (int i : filtered) filteredLines.push_back(splitlines[i]);
		//if (i < splitlines.size() - 4) filteredSides.push_back(sides[i]);
		return filteredLines;
	}

	static std::vector<SplitSide> getSplitLinesInNode(RaySpaceTree* rst, Node* leaf, bool storeCombi = false, 
														std::vector<Ray>& validCombis = std::vector<Ray>()) {

		std::vector<SplitSide> splitLines;

		rst->getSplittingLinesInNodeWithSide(leaf, splitLines);
		if (!rst->alldir) {
			std::vector<Ray> boxSides = rst->model->boundingCube.getCubeSideSquare(rst->maindir).quadLines;
			for (int i = 0; i < 4; i++) {
				Ray boxside = boxSides[i];
				boxside.index = splitLines.size();
				SplitSide split;
				split.ray = boxside;
				split.id = splitLines.size();
				splitLines.push_back(split);
			}
		}
		return filterSplittingLines(rst, leaf, splitLines, storeCombi, validCombis);
	}



	static std::vector<Ray> splitLineCombis(RaySpaceTree* rst, Node* leaf) {
		std::vector<Ray> validCombis;
		getSplitLinesInNode(rst, leaf, true, validCombis);
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