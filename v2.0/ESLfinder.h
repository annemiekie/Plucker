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
	std::vector<ESLCandidate> eslsNoVis;
	std::vector<ESLCandidate> esls;
	std::vector<Edge*> silhouetteEdges;
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
			for (auto& esl : esls) esl.ray.get3DfromPlucker();
			return true;
		}
		return false;
	};

private:

	bool cache = false;
	Cache<std::vector<Line4>>* combiCache;

	RaySpaceTree* rst;
	Primitive* prim;
	Node* node;

	ESLChecker eslChecker;

	bool printAll;
	bool print;

	bool storeAllESLs;

	std::vector<Vertex*> silhouetteVertices;
	std::vector<SplitSide> splitLines;

	std::vector<std::vector<int>> combi1E_store;
	std::vector<std::vector<int>> combi2E_store;
	std::vector<std::vector<int>> combi3E_store;
	CombiConfigurations combiS;

	bool splitVertexFlag = false;

	bool findEsl() {

		// Check if triangle edges are same as splitting lines and if yes, if it lies on the correct side
		if (!checkEdgeSplittingDuplicates()) return false;

		std::vector<std::vector<int>> cc;

		// Check without occlusion to see if prim can be excluded
		if (!checkCombi("SSV(T)", combiS.c2.size() * 3, 2, 0, 0, 2, combiS.c2, 0, cc, false) &&
			!checkCombi("SSST", combiS.c3.size() * 3, 3, 0, 0, 1, combiS.c3, 0, cc, false) &&
			!checkCombi("SSSS", combiS.c4.size(), 4, 0, 0, 0, combiS.c4, 0, cc, false)) return false;
		// maybe these can be used to check for degeneracies?


		// Check basic combis of extremal stabbing lines
		if (checkCombi("SSV(T)", combiS.c2.size() * 3, 2, 0, 0, 2, combiS.c2)) return true;
		if (checkCombi("SSST", combiS.c3.size() * 3, 3, 0, 0, 1, combiS.c3)) return true;

		checkHardCombis = true;

		// Find silhouette edges for primitive
		std::vector<Edge*> silhouetteEdgesFirst;
		std::vector<Edge*> silhouetteEdgesSecond;

		findSilhouettes(silhouetteEdgesFirst, silhouetteEdgesSecond);

		// Check all combis involving silhouette edges of some sort
		std::set<Vertex*, Vertex::cmp_ptr> silhEdgeVertices;

		if (silhouetteEdgesFirst.size() > 0)  if (checkSilhouetteCombis(silhouetteEdgesFirst, silhEdgeVertices)) return true;

		// Check basic (but large) combi of extremal stabbing lines
		if (checkCombi("SSSS", combiS.c4.size(), 4, 0, 0, 0, combiS.c4)) return true;

		// Check combis involving second tier silhouette edges
		if (silhouetteEdgesSecond.size() > 0) if (checkSilhouetteCombis(silhouetteEdgesSecond, silhEdgeVertices)) return true;

		return storeAllESLs && esls.size();// && !checkDegenerateESLs();
	}

	void findSilhouettes(std::vector<Edge*>& silhouetteEdgesFirst, std::vector<Edge*>& silhouetteEdgesSecond) {
		std::vector<Edge*> foundsilhouettes;
		if (rst->filledExact)
			getSilhouetteEdges(foundsilhouettes, node->primitiveSet);
		if (node->parent->filledExact)
			getSilhouetteEdges(foundsilhouettes, node->parent->primitiveSet);
		else
			getSilhouetteEdges(foundsilhouettes);

		for (Edge* e : foundsilhouettes) {
			bool found = false;
			for (Primitive* pr : e->triangles) {
				if (node->containsPrim(pr->id)) {
					if (checkEdgeInLeafCombis(e, combiS.c2) || checkEdgeInLeafCombis(e, combiS.c3)) {
						silhouetteEdgesFirst.push_back(e);
						found = true; break;
					}
				}
			}
			if (!found && !rst->filledExact) {
				if (checkEdgeInLeafCombis(e, combiS.c2) || checkEdgeInLeafCombis(e, combiS.c3))
					silhouetteEdgesSecond.push_back(e);
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

						splitVertexFlag = true;
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
					if (e->intersectsRay(esl.ray)) return true;
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
				if (vischeck) esls.push_back(esl);
				else eslsNoVis.push_back(esl);
				if (!storeAllESLs) return true;
			}
		}
		if (cache) combiCache->storeValueAtLastKey(cacheLines);
		return foundRay;
	}

	bool checkVeVt(std::vector<Vertex*>& silhVertices, int offset) {
		if (printAll) std::cout << "V(e)V(t) combi's: " << silhVertices.size() * 3 << std::endl;
		for (int i = offset; i < silhVertices.size(); i++) {// Vertex * vs : silhVertices) {
			Vertex* vs = silhVertices[i];
			for (Vertex* vt : prim->vertices) {
				Line4 ray(vs->pos, vt->pos);
				ESLCandidate esl = { ray };
				esl.silhouetteVertices = { vs };
				esl.triangleVertex = { vt };
				if (eslChecker.isExtremalStabbing(esl)) {
					if (print) std::cout << "V(e)V(t)" << std::endl;
					esls.push_back(esl);
					if (!storeAllESLs) true;
				}
			}
		}
		return false;
	};

	bool checkVeVe(std::vector<Vertex*>& silhVertices, int offset) {
		//std::vector<std::vector<int>> combi = Combinations::combi2(silhouettesize);

		if (printAll) std::cout << "V(e)V(e) combi's: " << silhVertices.size() * (silhVertices.size() - 1) << std::endl;
		for (int i = offset; i < silhVertices.size(); i++) {
			Vertex* vs = silhVertices[i];
			for (Vertex* vt : silhVertices) {
				if (vs->id == vt->id) continue;
				Line4 ray(vs->pos, vt->pos);
				ESLCandidate esl = { ray };
				esl.silhouetteVertices = { vs, vt };
				if (eslChecker.isExtremalStabbing(esl)) {
					if (print) std::cout << "V(e)V(e)" << std::endl;
					esls.push_back(esl);
					if (!storeAllESLs) true;
				}
			}
		}
		return false;
	};

	bool checkSilhouetteCombis(std::vector<Edge*>& silhouetteEdgesToAdd, std::set<Vertex*, Vertex::cmp_ptr>& silhVertices) {

		int offset_edge = silhouetteEdges.size();
		int offset_vertex = silhouetteVertices.size();

		std::set<Vertex*, Vertex::cmp_ptr> silhVertices1;
		for (Edge* toAdd : silhouetteEdgesToAdd) {
			bool add = true;
			for (Split& split : splitLines) {
				if (split.edge == toAdd) {
					add = false;
					break;
				}
			}
			if (add) silhouetteEdges.push_back(toAdd);
			for (Vertex* v : toAdd->vertices) {
				bool skip = false;
				for (Vertex* vp : prim->vertices) {
					if (vp == v) {
						skip = true;
						break;
					}
				}
				if (!skip) silhVertices1.insert(v);
			}
		}
		if (offset_edge == 0) {
			silhVertices = silhVertices1;
			for (Vertex* v : silhVertices) silhouetteVertices.push_back(v);
		}
		else {
			for (Vertex* v : silhVertices1) {
				if (silhVertices.find(v) == silhVertices.end()) silhouetteVertices.push_back(v);
			}
		}

		std::vector<std::vector<int>> e2;
		
		if (checkVeVt(silhouetteVertices, offset_vertex)) return true;

		if (checkCombi("SV(E)T", combiS.c1.size() * silhouetteEdges.size() * 2 * 3, 1, 2, 0, 1, combiS.c1, offset_vertex, e2)) return true;
		if (checkVeVe(silhouetteVertices, offset_vertex)) return true;

		std::vector<std::vector<int>> combi1E = Combinations::combi1(silhouetteEdges.size(), offset_edge);

		if (combi1E.size() > 0) {
			if (checkCombi("SV(E)E", combiS.c1.size() * combi1E.size() * silhouetteVertices.size(), 1, 2, 1, 0, combiS.c1, 0, combi1E)) return true;

			if (checkCombi("SEV(T)", combiS.c1.size() * combi1E.size() * 3, 1, 0, 1, 2, combiS.c1, 0, combi1E)) return true;

			if (checkCombi("V(E)ET", combi1E.size() * silhouetteVertices.size() * 3, 0, 2, 1, 1, e2, 0, combi1E)) return true;
		}
		if (silhouetteVertices.size() - offset_vertex > 0)
			if (checkCombi("SSV(E)", combiS.c2.size() * (silhouetteVertices.size() - offset_vertex), 2, 2, 0, 0, combiS.c2, offset_vertex, e2)) return true;

		if (offset_edge > 0) {
			if (combi1E_store.size() > 0) {
				if (checkCombi("SV(E)E", combiS.c1.size() * combi1E_store.size() * (silhouetteVertices.size() - offset_vertex), 1, 2, 1, 0, combiS.c1, offset_vertex, combi1E_store)) return true;
				if (checkCombi("V(E)ET", combi1E_store.size() * (silhouetteVertices.size() - offset_vertex) * 3, 0, 2, 1, 1, e2, offset_vertex, combi1E_store)) return true;
			}
			if (combi2E_store.size() > 0)
				if (checkCombi("V(E)EE", combi2E_store.size() * (silhouetteVertices.size() - offset_vertex), 0, 2, 2, 0, e2, offset_vertex, combi2E_store)) return true;
		}
		combi1E_store = combi1E;

		std::vector<std::vector<int>> combi2E;
		getEECombis(combi2E, offset_edge);
		for (std::vector<int>& combi : combi2E) combi2E_store.push_back(combi);

		if (combi2E.size() > 0) {
			if (checkCombi("EEV(T)", combi2E.size() * 3, 0, 0, 2, 2, e2, 0, combi2E)) return true;

			std::vector<std::vector<int>> combi3E;
			getEEECombis(combi2E_store, combi3E, offset_edge);
			for (std::vector<int>& combi : combi3E) combi3E_store.push_back(combi);

			if (combi3E.size() > 0) {
				if (checkCombi("EEET", combi3E.size() * 3, 0, 0, 3, 1, e2, 0, combi3E)) return true;
				if (checkCombi("SEEE", combiS.c1.size() * combi3E.size(), 1, 0, 3, 0, combiS.c1, 0, combi3E)) return true;

				std::vector<std::vector<int>> combi4E = Combinations::combiAddSelective(silhouetteEdges.size(), combi3E_store, offset_edge);
				if (combi4E.size() > 0) if (checkCombi("EEEE", combi4E.size() * 3, 0, 0, 4, 0, e2, 0, combi4E)) return true;
			}


			if (checkCombi("V(E)EE", combi2E.size() * silhouetteVertices.size(), 0, 2, 2, 0, e2, 0, combi2E)) return true;

			if (checkCombi("SEET", combiS.c1.size() * combi2E.size() * 3, 1, 0, 2, 1, combiS.c1, 0, combi2E)) return true;
			if (checkCombi("SSEE", combiS.c2.size() * combi2E.size(), 2, 0, 2, 0, combiS.c2, 0, combi2E)) return true;
		}

		if (combi1E.size() > 0) {
			if (checkCombi("SSET", combiS.c2.size() * combi1E.size() * 3, 2, 0, 1, 1, combiS.c2, 0, combi1E)) return true;
			if (checkCombi("SSSE", combiS.c3.size() * combi1E.size(), 3, 0, 1, 0, combiS.c3, 0, combi1E)) return true;
		}

		return false;
	}

	void getSilhouetteEdges(std::vector<Edge*>& silhouettes, std::set<int>& checktris = std::set<int>()) {
		robin_hood::unordered_map<uint64_t, Edge*> edgesToCheck;
		if (checktris.size() != 0) 	for (int t : checktris) for (Edge* e : rst->model->triangles[t]->edges) edgesToCheck[e->getKey()] = e;
		else edgesToCheck = rst->model->edges;

		for (auto& check_edge : edgesToCheck) {
			Edge* e = check_edge.second;
			// If edge is not in potential silhouettes
			//if (!rst->alldir && !e->silhouette) continue;

			// Don't use same triangle as silhouette edge.
			if (e->triangles[0] == prim) continue;
			if (e->triangles.size() == 2 && e->triangles[1] == prim) continue;

			// check if edge lies on correct side of primitive plane
			Vertex* edgeV1 = e->vertices[0];
			Vertex* edgeV2 = e->vertices[1];
			bool edgeIntersectPlane = false;
			if (!prim->vertexOnPositiveSidePlane(edgeV1) || !prim->vertexOnPositiveSidePlane(edgeV2)) {
				edgeIntersectPlane = true;
				if (!prim->vertexOnPositiveSidePlane(edgeV1) && !prim->vertexOnPositiveSidePlane(edgeV2)) continue;
			}

			// check if edge can be seen through window
			bool cntn = false;
			for (Primitive* p : e->triangles) {
				if (glm::dot(rst->window->normal, p->normal) < 0 && !rst->window->intersectsPlaneFromLines(p->getRayVector())) cntn = true;
				else cntn = false;
			}
			if (cntn) continue;

			bool sideCheck = false;
			bool intersect = false;
			std::vector<glm::dvec3> intersectpts;
			bool first = true;
			bool silhouetteFound = false;

			for (Vertex* primV : prim->vertices) {
				bool side = false;
				silhouetteFound = silhouetteFound || e->isSilhouetteForVertex(primV, side, rst->window->normal);

				if (!rst->alldir) {
					for (Vertex* edgeV : e->vertices) {
						Ray r = Ray(edgeV->pos, primV->pos);
						intersect = intersect || rst->window->inBounds(r);
						glm::vec3 windowIntersect = rst->window->rayIntersection(r);
						if (edgeIntersectPlane && rst->window->rayIntersectionDepth(r) > 0) windowIntersect *= -INFINITY;
						intersectpts.push_back(windowIntersect);
					}
				}
				if (!silhouetteFound) {
					if (first) { first = false;  sideCheck = side; }
					else if (side != sideCheck) silhouetteFound = true;
				}
				if (silhouetteFound && (rst->alldir || intersect)) break;
			}

			if (silhouetteFound) {
				if (!intersect) {
					AxisAlignedSquare window(rst->window);
					AxisAlignedSquare boundingSquare(intersectpts, rst->window);
					intersect = !window.noIntersection(boundingSquare);
				}
				if (intersect) silhouettes.push_back(e);
			}
		}
	}

	void getEEECombis(std::vector<std::vector<int>>& combi2Edges, std::vector<std::vector<int>>& combi3Edges, int offset = 0) {
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

	void getEECombis(std::vector<std::vector<int>>& combi2Edges, int offset = 0, int substract = 0) {
		if (silhouetteEdges.size() >= 2) {
			std::vector<std::vector<int>> combi = Combinations::combi2(silhouetteEdges.size(), offset, substract);

			for (std::vector<int>& c : combi) {
				Edge* e1 = silhouetteEdges[c[0]];
				Edge* e2 = silhouetteEdges[c[1]];
				if (rst->model->checkEdgeEdgeCache(e1, e2, rst->alldir, rst->maindir)) combi2Edges.push_back(c);
			}
		}
	}

	bool checkDegenerateESLs() {
		if (eslsNoVis.size() <= 2) return true;

		bool allintersect = true;
		for (int i = 0; i < eslsNoVis.size(); i++) {
			if (!eslsNoVis[i].ray.intersect(eslsNoVis[(i + 1) % eslsNoVis.size()].ray)) allintersect = false;
		}
		std::string tf = allintersect ? "true" : "false";
		if (node->index == 275 && prim->id == 246) std::cout << "intersect all? " << tf << std::endl;

		if (allintersect) return true;

		std::set<int> sameIds;
		for (Edge* e1 : eslsNoVis[0].allEdges) {
			for (Edge* e2 : eslsNoVis[1].allEdges) {
				if (e1->id == e2->id) sameIds.insert(e1->id);
			}
		}
		std::vector<int> notfound;
		for (int i = 2; i < eslsNoVis.size(); i++) {
			for (int id : sameIds) {
				bool found = false;

				for (Edge* e1 : eslsNoVis[i].allEdges) {
					if (e1->id == id) found = true;
				}
				if (!found) notfound.push_back(id);
			}
			for (int id : notfound) sameIds.erase(id);
		}
		tf = allintersect ? "true" : "false";
		if (node->index == 275 && prim->id == 246) std::cout << "same ids nr? " << sameIds.size() << std::endl;
		
		if (sameIds.size() == 0) return false;


		return true;
	}

	bool checkCombi(std::string combi_text, int combiNr, 
					int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges,
					std::vector<std::vector<int>>& splitLineCombis, int vertex_offset = 0,
					std::vector<std::vector<int>>& silhLineCombis = std::vector<std::vector<int>>(),
					bool vischeck = true) {

		if (printAll) std::cout << combi_text << " combi's: " << combiNr << std::endl;

		std::vector<Vertex*> vertexEdgeCheck;

		for (int g = 0; g < std::max(1, std::min(nrOfTriEdges, 1) * 3); g++) {
			if (nrOfTriEdges == 2) vertexEdgeCheck = { prim->vertices[(g + 1) % 3] };
			//if (nrOfTriEdges == 1) vertexEdgeCheck = { (int)model->indices[prim * 3 + edges[g][0]] , (int)model->indices[prim * 3 + edges[g][1]] };

			for (int h = vertex_offset; h < (nrOfVertices ? silhouetteVertices.size() : 1); h++) {
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
						//cntn = false;
						//for (int x = 0; x < esl.lines4.size(); x++) {
						//	for (int y = x + 1; y < esl.lines4.size(); y++) {
						//		if (esl.lines4[x].equal(esl.lines4[y])) cntn = true;
						//		// also check for shared vertices among the triedges/silhouette edges and splitlines -- don't check these!
						//		// need more info in eslcandidate
						//	}
						//}
						//if (cntn) continue;

						if (nrOfVertices) esl.silhouetteVertices = vertexEdgeCheck;
						if (nrOfTriEdges == 2) esl.triangleVertex.push_back(prim->vertices[(g + 1) % 3]);
						esl.fillUpLines();
						if (esl.combinationDegenerate()) continue;

						for (SplitSide& s : esl.splittingLines) {
								if (s.edge != NULL) {
									for (Edge* e : silhouetteEdges) {
										if (e->id == s.edge->id) s.edgeIsSilhouette = true;
									}
							}
						}

						if (checkRaysThrough4Lines(esl, vischeck)) { 
							if (print) {
								std::cout << combi_text << " ";
								if (!vischeck) std::cout << "No visibility Check: ";
								for (Ray& r : esl.lines4) std::cout << r.index << " ";
								std::cout << std::endl;
							}
							
							if (!storeAllESLs) {
								if (vischeck) return true;
								else if (!checkDegenerateESLs()) return true;
							}
							else if (!vischeck && !checkDegenerateESLs()) return true;
						}
					}
				}
			}
		}
		return false;
	};




};