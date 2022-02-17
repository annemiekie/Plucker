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
	Node* leaf;
	ESLChecker eslChecker;

	bool printAll;
	bool print;

	bool allESLs;
	std::vector<Line4> allESLsForPrimNode;

	std::vector<Edge*> silhouetteEdges;
	std::vector<Vertex*> silhouetteVertices;
	std::vector<Split> splitLines;

	std::vector<std::vector<int>> combi2E_store;
	CombiConfigurations combiS;
	bool checkAllPrimsForSilh = true;

	Line4 ray;

	ESLFinder() {};

	ESLFinder(RaySpaceTree* rst, Primitive* prim, Node* leaf,
				bool detail, bool print, bool allESLs, 
				std::vector<Split>& splitLines, CombiConfigurations& splitcombi,
				bool checkAllPrimsForSilh,
				bool cache = false, Cache<std::vector<Line4>>* combiCache = NULL) :

				rst(rst), prim(prim), leaf(leaf), print(print), printAll(detail), allESLs(allESLs),
				splitLines(splitLines), combiS(splitcombi),
				checkAllPrimsForSilh(checkAllPrimsForSilh), cache(cache), combiCache(combiCache)
	{	
		eslChecker = ESLChecker(leaf, prim, rst, printAll);
	};

	bool find() {

		bool printAll = false;
		int splitLinesSize = splitLines.size();

		// Check if triangle edges are same as splitting lines and if yes, if it lies on the correct side
		if (!checkEdgeSplittingDuplicates()) return false;

		std::vector<std::vector<int>> cc;

		// Check without occlusion to see if prim can be excluded
		if (!checkCombi("SSV(T)", combiS.c2.size() * 3, 2, 0, 0, 2, combiS.c2, cc, false) &&
			!checkCombi("SSST", combiS.c3.size() * 3, 3, 0, 0, 1, combiS.c3, cc, false) &&
			!checkCombi("SSSS", combiS.c4.size(), 4, 0, 0, 0, combiS.c4, cc, false)) return false;

		// Check basic combis of extremal stabbing lines
		if (checkCombiStoreAll("SSV(T)", combiS.c2.size() * 3, 2, 0, 0, 2, combiS.c2)) return true;
		if (checkCombiStoreAll("SSST", combiS.c3.size() * 3, 3, 0, 0, 1, combiS.c3)) return true;

		// Find silhouette edges for primitive
		std::vector<Edge*> silhouetteEdges;
		std::vector<Edge*> silhouetteEdgesFirst;
		std::vector<Edge*> silhouetteEdgesSecond;

		if (checkAllPrimsForSilh) {
			rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, silhouetteEdges);
			//exact would be to test if triangles belong in leaf
			for (Edge* e : silhouetteEdges) {
				bool found = false;
				for (Primitive* pr : e->triangles) {
					if (leaf->primitiveSet.find(pr->id) != leaf->primitiveSet.end()) {
						silhouetteEdgesFirst.push_back(e);
						//if (print) std::cout << "Silh Edge " << e->id << std::endl;
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
		else rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, silhouetteEdgesFirst, leaf->parent->primitiveSet);

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
		silhouetteEdges = std::vector<Edge*>();

		if (silhouetteEdgesFirst.size() > 0) {
			if (checkSilhouetteCombis(silhouetteEdgesFirst, silhEdgeVertices)) return true;
		}

		// Check basic (but large) combi of extremal stabbing lines
		if (checkCombiStoreAll("SSSS", combiS.c4.size(), 4, 0, 0, 0, combiS.c4)) return true;

		// Check combis involving second tier silhouette edges
		if (silhouetteEdgesSecond.size() > 0) if (checkSilhouetteCombis(silhouetteEdgesSecond, silhEdgeVertices)) return true;

		if (allESLs && allESLsForPrimNode.size()) return true;
		return false;
	}


	bool checkEdgeSplittingDuplicates() {
		for (SplitSide& split : splitLines) {
			bool throughEdge = false;
			for (Edge* e : prim->edges) {
				if (split.line.equal(e->ray, 1E-5)) { // not very precise?!
					throughEdge = true;
					if (!prim->onRightSideOfSplitEdge(split.line, split.side)) return false;
				}
			}
			if (!throughEdge) {
				for (Vertex* v : prim->vertices) {
					if (split.line.throughVertex(v)) {
						if (!prim->onRightSideOfSplitVertex(v, split.line, split.side)) return false;
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
				for (int i = 0; i < combi[0].size(); i++) lines4.push_back(splitLines[combi[c][i]].line);
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
	bool checkRaysThrough4Lines(std::vector<Ray>& lines4, bool vischeck, 
								std::vector<uint64_t>& indices = std::vector<uint64_t>(),
								std::vector<Vertex*>& silhVertices = std::vector<Vertex*>(), 
								std::vector<Edge*>& silhEdges = std::vector<Edge*>()) {

		std::vector<Line4> intersectLines;
		if (!cache || !combiCache->getValue(indices, intersectLines)) intersectLines = Lines4Finder::find(lines4);
		std::vector<Line4> cacheLines;

		for (int i = 0; i < intersectLines.size(); i++) {
			ray = intersectLines[i];
			bool inBox = true;
			ESLCandidate esl = { ray, lines4, silhVertices, silhEdges };
			bool checkRay = eslChecker.isExtremalStabbing(esl, vischeck);
			if (cache && esl.inBox) {
				ray.checkboth = false;
				cacheLines.push_back(ray);
			}
			if (checkRay) {
				if (cache) {
					if (i == 0) cacheLines.push_back(intersectLines[1]);
					combiCache->storeValueAtLastKey(cacheLines);
				}
				return true;
			}
		}
		if (cache) combiCache->storeValueAtLastKey(cacheLines);
		return false;
	}

	bool checkVeVt(std::vector<Vertex*>& silhVertices) {
		if (printAll) std::cout << "V(e)V(t) combi's: " << silhVertices.size() * 3 << std::endl;
		for (Vertex* vs : silhVertices) {
			for (Vertex* vt : prim->vertices) {
				ray = Line4(vs->pos, vt->pos);
				ESLCandidate esl = { ray, std::vector<Ray>(), { vs } };
				if (eslChecker.isExtremalStabbing(esl)) {
					if (print) std::cout << "V(e)V(t)" << std::endl;
					if (allESLs) allESLsForPrimNode.push_back(ray);
					else return true;
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
				ray = Line4(vs->pos, vt->pos);
				ESLCandidate esl = { ray, std::vector<Ray>(), { vs, vt } };
				if (eslChecker.isExtremalStabbing(esl)) {
					if (print) std::cout << "V(e)V(e)" << std::endl;
					if (allESLs) allESLsForPrimNode.push_back(ray);
					else return true;
				}
			}
		}
		return false;
	};

	bool checkCombiStoreAll(std::string combi_text, int combiNr, int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges,
							std::vector<std::vector<int>>& splitLineCombis, std::vector<std::vector<int>>& silhLineCombis = std::vector<std::vector<int>>(),
							bool vischeck = true) {
		return checkCombi(combi_text, combiNr, nrOfsplitLines, nrOfVertices, nrOfsilhEdges, nrOfTriEdges,
							splitLineCombis, silhLineCombis, vischeck) && !allESLs;
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

		std::vector<Vertex*> silhouetteVertices;
		for (Vertex* v : silhVertices) silhouetteVertices.push_back(v);

		std::vector<std::vector<int>> e2;

		if (checkVeVt(silhouetteVertices) && !allESLs) return true;

		if (checkCombiStoreAll("SV(E)T", combiS.c1.size() * silhouetteEdges.size() * 2 * 3, 1, 2, 0, 1, combiS.c1, e2)) return true;
		if (checkVeVe(silhouetteVertices) && !allESLs) return true;

		std::vector<std::vector<int>> combi1E = Combinations::combi1(silhouetteEdges.size(), offset);

		if (checkCombiStoreAll("SV(E)E", combiS.c1.size() * combi1E.size() * silhouetteEdges.size() * 2, 1, 2, 1, 0, combiS.c1, combi1E)) return true;
		if (checkCombiStoreAll("SEV(T)", combiS.c1.size() * combi1E.size() * 3, 1, 0, 1, 2, combiS.c1, combi1E)) return true;

		if (checkCombiStoreAll("SSV(E)", combiS.c2.size() * silhouetteEdges.size() * 2, 2, 2, 0, 0, combiS.c2, e2)) return true;
		if (checkCombiStoreAll("V(E)ET", combi1E.size() * silhouetteEdges.size() * 2 * 3, 0, 2, 1, 1, e2, combi1E)) return true;

		std::vector<std::vector<int>> combi2E;
		getEECombis(silhouetteEdges, combi2E, offset);
		for (std::vector<int>& combi : combi2E) combi2E_store.push_back(combi);

		if (combi2E.size() > 0) {
			if (checkCombiStoreAll("EEV(T)", combi2E.size() * 3, 0, 0, 2, 2, e2, combi2E)) return true;

			std::vector<std::vector<int>> combi3E;
			getEEECombis(silhouetteEdges, combi2E_store, combi3E, offset);

			if (combi3E.size() > 0) {
				if (checkCombiStoreAll("EEET", combi3E.size() * 3, 0, 0, 3, 1, e2, combi3E)) return true;
				if (checkCombiStoreAll("SEEE", combiS.c1.size() * combi3E.size(), 1, 0, 3, 0, combiS.c1, combi3E)) return true;

				std::vector<std::vector<int>> combi4E = Combinations::combiAddSelective(silhouetteEdges.size(), combi3E, offset);
				if (combi4E.size() > 0) {
					if (checkCombiStoreAll("EEEE", combi4E.size() * 3, 0, 0, 4, 0, e2, combi4E)) return true;
				}
			}

			if (checkCombiStoreAll("V(E)EE", combi2E.size() * silhouetteEdges.size() * 2, 0, 2, 2, 0, e2, combi2E)) return true;
			if (checkCombiStoreAll("SEET", combiS.c1.size() * combi2E.size() * 3, 1, 0, 2, 1, combiS.c1, combi2E)) return true;
			if (checkCombiStoreAll("SSEE", combiS.c2.size() * combi2E.size(), 2, 0, 2, 0, combiS.c2, combi2E)) return true;
		}

		if (checkCombiStoreAll("SSET", combiS.c2.size() * combi1E.size() * 3, 2, 0, 1, 1, combiS.c2, combi1E)) return true;
		if (checkCombiStoreAll("SSSE", combiS.c3.size() * combi1E.size(), 3, 0, 1, 0, combiS.c3, combi1E)) return true;

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

			for (int h = 0; h < std::max(1, (int)silhouetteVertices.size()); h++) {
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

						std::vector<Ray> lines4;
						std::vector<Edge*> silhEdges;
						std::vector<Vertex*> silhVertices;
						std::vector<uint64_t> indices;

						int linecount = 0;
						for (int k = 0; k < nrOfsplitLines; k++) {
							lines4.push_back(splitLines[splitLineCombis[j][k]].line);
							//	if (cacheCombi) indices.push_back(lines[linecount].index + model->edges.size() + model->vertices.size());
							linecount++;
						}

						for (int k = 0; k < nrOfsilhEdges; k++) {
							lines4.push_back(silhouetteEdges[silhLineCombis[i][k]]->ray);
							silhEdges.push_back(silhouetteEdges[silhLineCombis[i][k]]);
							//	if (cacheCombi) indices.push_back(lines[linecount].index);
							linecount++;
						}


						for (int k = 0; k < nrOfTriEdges; k++) {
							lines4.push_back(prim->rays[(g + k) % 3]);// triEdgeRays[edges[g][k]]);
						//	if (cacheCombi) indices.push_back(lines[linecount].index);
							linecount++;
						}

						// check equal lines
						for (int x = 0; x < linecount; x++) for (int y = x + 1; y < linecount; y++)
							if (lines4[x].equal(lines4[y])) continue;

						cntn = false;
						if (nrOfVertices) {
							for (Ray& r : lines4) if (r.throughVertex(silhouetteVertices[h])) cntn = true;
							if (cntn) continue;
							lines4.push_back(silhouetteVertices[h]->edges[0]->ray);
							lines4.push_back(silhouetteVertices[h]->edges[1]->ray);
							silhVertices = vertexEdgeCheck;
							//	if (cacheCombi) {
							//		indices.push_back(lines[linecount].index + model->edges.size());
							//		indices.push_back(lines[linecount + 1].index + model->edges.size());
							//	}
							linecount += 2;
						}

						if (checkRaysThrough4Lines(lines4, vischeck, indices, silhVertices, silhEdges)) { // change this!!
						//if (checkRaysThroughLines(rst, prim, ray, leaf, printAll, lines4, indices, vischeck, vertexEdgeCheck, silhEdges)) { // change this!!
							if (print) {
								std::cout << combi_text << " ";
								for (Ray& r : lines4) std::cout << r.index << " ";
								std::cout << std::endl;
							}
							if (allESLs) allESLsForPrimNode.push_back(ray);
							else return true;
						}
					}
				}
			}
		}
		return false;
	};




};