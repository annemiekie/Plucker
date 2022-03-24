#pragma once
#include "cache.h"
#include "model.h"
#include "node.h"
#include "findLineThroughFour.h"
#include "combinations.h"
#include "rst.h"
#include "eslChecker.h"
#include "eslCandidate.h"
#include "splitSide.h"
#include "combiConfigurations.h"

class ESLFinder {
public:
	bool cache = false;
	Cache<std::vector<Line4>>* combiCache;

	RaySpaceTree* rst;
	Primitive* prim;
	Node* node;

	ESLChecker eslChecker;

	bool printAll;
	bool print;

	bool storeAllESLs;
	std::vector<Line4> esls;

	std::vector<Edge*> silhouetteEdges;
	std::vector<Vertex*> silhouetteVertices;
	std::vector<SplitSide> splitLines;

	std::vector<std::vector<int>> combi2E_store;
	std::vector<std::vector<int>> combi3E_store;
	CombiConfigurations combiS;

	bool checkHardCombis = false;

	ESLFinder() {};

	ESLFinder(RaySpaceTree* rst, Primitive* prim, Node* node,
				bool detail, bool print, bool allESLs, 
				std::vector<SplitSide>& splitLines, CombiConfigurations& splitcombi,
				bool cache = false, Cache<std::vector<Line4>>* combiCache = NULL) :

				rst(rst), prim(prim), node(node), print(print), printAll(detail), storeAllESLs(allESLs),
				splitLines(splitLines), combiS(splitcombi), cache(cache), combiCache(combiCache)
	{	
		eslChecker = ESLChecker(node, prim, rst, printAll);
	};

	bool find() {
		if (findEsl()) {
			for (auto& r : esls) r.get3DfromPlucker();
			return true;
		}
		return false;
	};

	bool findEsl() {

		// Check if triangle edges are same as splitting lines and if yes, if it lies on the correct side
		if (!checkEdgeSplittingDuplicates()) return false;

		std::vector<std::vector<int>> cc;

		// Check without occlusion to see if prim can be excluded
		if (!checkCombi("SSV(T)", combiS.c2.size() * 3, 2, 0, 0, 2, combiS.c2, cc, false) &&
			!checkCombi("SSST", combiS.c3.size() * 3, 3, 0, 0, 1, combiS.c3, cc, false) &&
			!checkCombi("SSSS", combiS.c4.size(), 4, 0, 0, 0, combiS.c4, cc, false)) return false;

		// Check basic combis of extremal stabbing lines
		if (checkCombi("SSV(T)", combiS.c2.size() * 3, 2, 0, 0, 2, combiS.c2)) return true;
		if (checkCombi("SSST", combiS.c3.size() * 3, 3, 0, 0, 1, combiS.c3)) return true;

		checkHardCombis = true;

		// Find silhouette edges for primitive
		std::vector<Edge*> silhouetteEdgesFirst;
		std::vector<Edge*> silhouetteEdgesSecond;

		findSilhouettes(silhouetteEdgesFirst, silhouetteEdgesSecond);

		//if (getedges) {
		//	for (Edge* e : silhouetteEdgesFirst) {
		//		for (Vertex* v : e->vertices) edges.push_back(v->pos);
		//	}
		//	for (Edge* e : silhouetteEdgesSecond) {
		//		for (Vertex* v : e->vertices) edges.push_back(v->pos);
		//	}
		//}

		// Check all combis involving silhouette edges of some sort
		std::set<Vertex*, Vertex::cmp_ptr> silhEdgeVertices;

		if (silhouetteEdgesFirst.size() > 0) {
			if (checkSilhouetteCombis(silhouetteEdgesFirst, silhEdgeVertices)) return true;
		}

		// Check basic (but large) combi of extremal stabbing lines
		if (checkCombi("SSSS", combiS.c4.size(), 4, 0, 0, 0, combiS.c4)) return true;

		// Check combis involving second tier silhouette edges
		if (silhouetteEdgesSecond.size() > 0) if (checkSilhouetteCombis(silhouetteEdgesSecond, silhEdgeVertices)) return true;

		return esls.size();
	}

	void findSilhouettes(std::vector<Edge*>& silhouetteEdgesFirst, std::vector<Edge*>& silhouetteEdgesSecond) {
		std::vector<Edge*> foundsilhouettes;
		if (rst->filledExact) {
			rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, silhouetteEdgesFirst, node->primitiveSet);
			return;
		}
		else if (node->parent->filledExact)
			rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, foundsilhouettes, node->parent->primitiveSet);
		else
			rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, foundsilhouettes);

		for (Edge* e : foundsilhouettes) {
			bool found = false;
			for (Primitive* pr : e->triangles) {
				if (node->primitiveSet.find(pr->id) != node->primitiveSet.end()) {
					silhouetteEdgesFirst.push_back(e);
					found = true; break;
				}
			}
			if (!found) {
				Line4 eslEdge;
				if (checkEdgeInLeafCombis(e, combiS.c2) || checkEdgeInLeafCombis(e, combiS.c3)) {
					silhouetteEdgesSecond.push_back(e);
					//if (getedges) {
					//	eslEdge.get3DfromPlucker();
					//	eslEdges.push_back(eslEdge);
					//}
				}
			}
		}
	}

	bool checkEdgeSplittingDuplicates() {
		for (SplitSide& split : splitLines) {
			bool throughEdge = false;
			for (Edge* e : prim->edges) {
				if (split.ray.equal(e->ray, 1E-5)) { // not very precise?!
					throughEdge = true;
					if (!prim->onRightSideOfSplitEdge(split.ray, split.side)) return false;
				}
			}
			if (!throughEdge) {
				for (Vertex* v : prim->vertices) {
					if (split.ray.throughVertex(v)) {
						if (prim->onRightSideOfSplitVertex(v, split.ray, split.side)) break;

						Vertex* otherVertex;
						// Check which of the edges attached to the vertex is the splitline
						for (Edge* e : v->edges) {
							if (e->ray.equal(split.ray, 1E-5)) {
								for (Vertex* v2 : e->vertices) if (v != v2) otherVertex = v2;
							}
						}

						if (prim->getPlane().pointOnPlane(otherVertex->pos)) return false;
						else if (prim->getPlane().pointOnPositiveSide(otherVertex->pos)) return false;
						break;
					}
				}
			}
		}
		return true;
	}

	bool checkEdgeInLeafCombis(Edge* e, std::vector<std::vector<int>>& combi) {

		std::vector<std::vector<Ray>> edgeCombi;
		if (combi[0].size() == 2) edgeCombi = { {e->vertices[0]->edges[0]->ray, e->vertices[0]->edges[1]->ray},
												{e->vertices[1]->edges[0]->ray, e->vertices[1]->edges[1]->ray} };
		else if (combi[0].size() == 3) edgeCombi = { {e->ray} };
		else return false;
		
		for (std::vector<Ray>& rays : edgeCombi) {
			for (int c = 0; c < combi.size(); c++) {
				std::vector<Ray> lines4;
				for (int i = 0; i < combi[0].size(); i++) lines4.push_back(splitLines[combi[c][i]].ray);
				for (Ray& r : rays) lines4.push_back(r);
		
				std::vector<Line4> intersectLines = Lines4Finder::find(lines4);		
				for (Line4& r: intersectLines) {
					bool inbox = false;
					ESLCandidate esl = { r, lines4 };
					if (!eslChecker.isExtremalStabbing(esl, false, false)) continue;
					if (combi[0].size() == 2) return true;
					float t = e->ray.depthToIntersectionWithRay(r);
					if (t > 0 && t < 1) return true;
				}
			}
		}
		return false;
	}

	// should store potentially two lines?? for all ESLs
	bool checkRaysThrough4Lines(ESLCandidate& esl, bool vischeck) {

		std::vector<Line4> intersectLines;
		if (!cache || !combiCache->getValue(esl.indices, intersectLines)) intersectLines = Lines4Finder::find(esl.lines4);
		std::vector<Line4> cacheLines;
		bool foundRay = false;

		for (int i = 0; i < intersectLines.size(); i++) {
			bool inBox = true;
			esl.ray = intersectLines[i];
			bool checkRay = eslChecker.isExtremalStabbing(esl, vischeck);
			if (cache && esl.inBox) {
				esl.ray.checkboth = false;
				cacheLines.push_back(esl.ray);
			}
			if (checkRay) {
				if (cache) {
					if (i == 0) cacheLines.push_back(intersectLines[1]);
					combiCache->storeValueAtLastKey(cacheLines);
				}
				foundRay = true;
				if (vischeck) esls.push_back(esl.ray);
				if (!storeAllESLs) return true;
			}
		}
		if (cache) combiCache->storeValueAtLastKey(cacheLines);
		return foundRay;
	}

	bool checkVeVt(std::vector<Vertex*>& silhVertices) {
		if (printAll) std::cout << "V(e)V(t) combi's: " << silhVertices.size() * 3 << std::endl;
		for (Vertex* vs : silhVertices) {
			for (Vertex* vt : prim->vertices) {
				Line4 ray(vs->pos, vt->pos);
				ESLCandidate esl = { ray };
				esl.silhouetteVertices = { vs };
				if (eslChecker.isExtremalStabbing(esl)) {
					if (print) std::cout << "V(e)V(t)" << std::endl;
					esls.push_back(ray);
					if (!storeAllESLs) true;
				}
			}
		}
		return false;
	};

	bool checkVeVe(std::vector<Vertex*>& silhVertices) {
		//std::vector<std::vector<int>> combi = Combinations::combi2(silhouettesize);

		if (printAll) std::cout << "V(e)V(e) combi's: " << silhVertices.size() * (silhVertices.size() - 1) << std::endl;
		for (Vertex* vs : silhVertices) {
			for (Vertex* vt : silhVertices) {
				if (vs->id == vt->id) continue;
				Line4 ray(vs->pos, vt->pos);
				ESLCandidate esl = { ray };
				esl.silhouetteVertices = { vs, vt };
				if (eslChecker.isExtremalStabbing(esl)) {
					if (print) std::cout << "V(e)V(e)" << std::endl;
					esls.push_back(ray);
					if (!storeAllESLs) true;
				}
			}
		}
		return false;
	};

	bool checkCombiStoreAll(std::string combi_text, int combiNr, int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges,
							std::vector<std::vector<int>>& splitLineCombis, std::vector<std::vector<int>>& silhLineCombis = std::vector<std::vector<int>>(),
							bool vischeck = true) {
		return checkCombi(combi_text, combiNr, nrOfsplitLines, nrOfVertices, nrOfsilhEdges, nrOfTriEdges,
							splitLineCombis, silhLineCombis, vischeck) && !storeAllESLs;
	}

	bool checkSilhouetteCombis(std::vector<Edge*>& silhouetteEdgesToAdd, std::set<Vertex*, Vertex::cmp_ptr>& silhVertices) {


		for (Edge* toAdd : silhouetteEdgesToAdd) {
			silhouetteEdges.push_back(toAdd);
			for (Vertex* v : toAdd->vertices) {
				bool skip = false;
				for (Vertex* vp : prim->vertices) {
					if (vp == v) {
						skip = true;
						break;
					}
				}
				if (!skip) silhVertices.insert(v);
			}
		}

		int offset = silhouetteEdges.size() - silhouetteEdgesToAdd.size();

		for (Vertex* v : silhVertices) silhouetteVertices.push_back(v);

		std::vector<std::vector<int>> e2;

		if (checkVeVt(silhouetteVertices)) return true;

		if (checkCombi("SV(E)T", combiS.c1.size() * silhouetteEdges.size() * 2 * 3, 1, 2, 0, 1, combiS.c1, e2)) return true;
		if (checkVeVe(silhouetteVertices)) return true;

		std::vector<std::vector<int>> combi1E = Combinations::combi1(silhouetteEdges.size(), offset);

		if (checkCombi("SV(E)E", combiS.c1.size() * combi1E.size() * silhouetteEdges.size() * 2, 1, 2, 1, 0, combiS.c1, combi1E)) return true;
		if (checkCombi("SEV(T)", combiS.c1.size() * combi1E.size() * 3, 1, 0, 1, 2, combiS.c1, combi1E)) return true;

		if (checkCombi("SSV(E)", combiS.c2.size() * silhouetteEdges.size() * 2, 2, 2, 0, 0, combiS.c2, e2)) return true;
		if (checkCombi("V(E)ET", combi1E.size() * silhouetteEdges.size() * 2 * 3, 0, 2, 1, 1, e2, combi1E)) return true;

		std::vector<std::vector<int>> combi2E;
		getEECombis(silhouetteEdges, combi2E, offset);
		for (std::vector<int>& combi : combi2E) combi2E_store.push_back(combi);

		if (combi2E.size() > 0) {
			if (checkCombi("EEV(T)", combi2E.size() * 3, 0, 0, 2, 2, e2, combi2E)) return true;

			std::vector<std::vector<int>> combi3E;
			getEEECombis(silhouetteEdges, combi2E_store, combi3E, offset);
			for (std::vector<int>& combi : combi3E) combi3E_store.push_back(combi);

			if (combi3E.size() > 0) {
				if (checkCombi("EEET", combi3E.size() * 3, 0, 0, 3, 1, e2, combi3E)) return true;
				if (checkCombi("SEEE", combiS.c1.size() * combi3E.size(), 1, 0, 3, 0, combiS.c1, combi3E)) return true;

				std::vector<std::vector<int>> combi4E = Combinations::combiAddSelective(silhouetteEdges.size(), combi3E_store, offset);
				if (combi4E.size() > 0) if (checkCombi("EEEE", combi4E.size() * 3, 0, 0, 4, 0, e2, combi4E)) return true;
			}

			if (checkCombi("V(E)EE", combi2E.size() * silhouetteEdges.size() * 2, 0, 2, 2, 0, e2, combi2E)) return true;
			if (checkCombi("SEET", combiS.c1.size() * combi2E.size() * 3, 1, 0, 2, 1, combiS.c1, combi2E)) return true;
			if (checkCombi("SSEE", combiS.c2.size() * combi2E.size(), 2, 0, 2, 0, combiS.c2, combi2E)) return true;
		}

		if (checkCombi("SSET", combiS.c2.size() * combi1E.size() * 3, 2, 0, 1, 1, combiS.c2, combi1E)) return true;
		if (checkCombi("SSSE", combiS.c3.size() * combi1E.size(), 3, 0, 1, 0, combiS.c3, combi1E)) return true;

		return false;
	}


	void getEEECombis(std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges, 
						std::vector<std::vector<int>>& combi3Edges, int offset) {
		if (silhouetteEdges.size() >= 3 && combi2Edges.size() > 0) {
			std::vector<std::vector<int>> combi = Combinations::combiAddSelective(silhouetteEdges.size(), combi2Edges, offset);
			for (std::vector<int>& c : combi) {
				Edge* e1 = silhouetteEdges[c[0]];
				Edge* e2 = silhouetteEdges[c[1]];
				Edge* e3 = silhouetteEdges[c[2]];

				if (rst->model->checkEdgeEdgeEdge(e1, e2, e3)) combi3Edges.push_back(c);
			}
		}
	}

	void getEECombis(std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges, int offset) {
		if (silhouetteEdges.size() >= 2) {
			std::vector<std::vector<int>> combi = Combinations::combi2(silhouetteEdges.size(), offset);

			for (std::vector<int>& c : combi) {
				Edge* e1 = silhouetteEdges[c[0]];
				Edge* e2 = silhouetteEdges[c[1]];
				if (rst->model->checkEdgeEdgeCache(e1, e2, rst->alldir, rst->maindir)) combi2Edges.push_back(c);
			}
		}
	}


	bool checkCombi(std::string combi_text, int combiNr, 
					int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges,
					std::vector<std::vector<int>>& splitLineCombis, 
					std::vector<std::vector<int>>& silhLineCombis = std::vector<std::vector<int>>(),
					bool vischeck = true) {

		if (printAll) std::cout << combi_text << " combi's: " << combiNr << std::endl;

		std::vector<Vertex*> vertexEdgeCheck;

		for (int g = 0; g < std::max(1, std::min(nrOfTriEdges, 1) * 3); g++) {
			if (nrOfTriEdges == 2) vertexEdgeCheck = { prim->vertices[(g + 1) % 3] };
			//if (nrOfTriEdges == 1) vertexEdgeCheck = { (int)model->indices[prim * 3 + edges[g][0]] , (int)model->indices[prim * 3 + edges[g][1]] };

			for (int h = 0; h < (nrOfVertices ? silhouetteVertices.size() : 1); h++) {
				if (nrOfVertices > 0) vertexEdgeCheck = { silhouetteVertices[h] };

				for (int i = 0; i < std::max(1, (int)silhLineCombis.size()); i++) {
					bool cntn = false;
					for (int ix = 0; ix < nrOfsilhEdges; ix++) {
						cntn = true;
						if (nrOfsilhEdges > 0 && vertexEdgeCheck.size() == 1) {
							Edge* e = silhouetteEdges[silhLineCombis[i][ix]];
							// if silhouette edge vertices are the same as one of the vertices used break
							if (vertexEdgeCheck[0] == e->vertices[0] || vertexEdgeCheck[0] == e->vertices[1]) break;
							//glm::vec3 vert = rst->model->verticesIndexed[vertexEdgeCheck[0]];
							bool side;
							if (!e->isSilhouetteForPos(vertexEdgeCheck[0]->pos, side)) break;
							// rst->model->checkSilhouetteEdge(vert, e, true, glm::vec3(0), side) <= 0) break;
							//extra check
							//if (!model->boundingCube.intersectSideSwath(vert, model->verticesIndexed[e.v[0]], model->verticesIndexed[e.v[1]], maindir)) break;
						}
						//if ((nrOfsilhEdges > 2 || (nrOfsilhEdges == 2 && ix == 0)) && vertexEdgeCheck.size() > 0) {
						//	std::vector<glm::vec3> n;
						//	std::vector<float> d;
						//	spaceSpannedByEdges(silhouetteEdges[silhLineCombis[i][ix]], silhouetteEdges[silhLineCombis[i][(ix + 1) % nrOfsilhEdges]], n, d);
						//	std::vector<glm::vec3> checkPoints;
						//	for (int vec : vertexEdgeCheck) checkPoints.push_back(model->verticesIndexed[vec]);
						//	if (!checkPointsInHalfSpaces(n, d, checkPoints)) break;
						//}
						cntn = false;
					}
					if (cntn) continue;

					for (int j = 0; j < std::max(1, (int)splitLineCombis.size()); j++) {

						ESLCandidate esl;
						//std::vector<uint64_t> indices;

						int linecount = 0;
						for (int k = 0; k < nrOfsplitLines; k++) esl.splittingLines.push_back(splitLines[splitLineCombis[j][k]]);
						for (int k = 0; k < nrOfsilhEdges; k++) esl.silhouetteEdges.push_back(silhouetteEdges[silhLineCombis[i][k]]);
						for (int k = 0; k < nrOfTriEdges; k++) esl.triangleEdges.push_back(prim->edges[(g + k) % 3]);

						// check equal lines
						cntn = false;
						for (int x = 0; x < esl.lines4.size(); x++)
							for (int y = x + 1; y < esl.lines4.size(); y++) 
								if (esl.lines4[x].equal(esl.lines4[y])) cntn = true;
						if (cntn) continue;

						if (nrOfVertices) {
							for (Ray& r : esl.lines4) if (r.throughVertex(silhouetteVertices[h])) cntn = true;
							if (cntn) continue;
							esl.silhouetteVertices = vertexEdgeCheck;
							//	if (cacheCombi) {
							//		indices.push_back(lines[linecount].index + model->edges.size());
							//		indices.push_back(lines[linecount + 1].index + model->edges.size());
							//	}
						}
						esl.fillUpLines();

						if (checkRaysThrough4Lines(esl, vischeck)) { 
							if (print && vischeck) {
								std::cout << combi_text << " ";
								for (Ray& r : esl.lines4) std::cout << r.index << " ";
								std::cout << std::endl;
							}
							if (!storeAllESLs || !vischeck) return true;
						}
					}
				}
			}
		}
		return false;
	};




};