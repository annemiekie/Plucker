#pragma once
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

class ESLFindAll {
public:
	std::vector<ESLCandidate> eslsNoVis;
	std::vector<ESLCandidate> esls;
	std::vector<Edge*> silhouetteEdges;
	bool checkHardCombis = false;

	ESLFindAll() {};

	ESLFindAll(RaySpaceTree* rst) :rst(rst)
	{
		//eslChecker = ESLChecker(node, prim, rst, printAll);
	};

	void find() {
		findAllEsl();
	};

private:

	RaySpaceTree* rst;

	ESLChecker eslChecker;

	std::vector<Vertex*> silhouetteVertices;
	std::vector<SplitSide> splitLines;

	//std::vector<std::vector<int>> E; // all possible silhouette edges for window AND splitline generators
	//std::vector<std::vector<int>> EEE;
	//std::vector<std::vector<int>> EEEE;

	std::set<int> T; // all viable edges (can be seen from window, depends on planes of primitives attached)
	std::set<int> Vt;
	std::set<std::vector<int>> ET;
	std::set<std::vector<int>> EE;

	std::set<std::vector<int>> VeT;
	std::set<std::vector<int>> EVe;
	std::set<std::vector<int>> EVt;

	std::set<std::vector<int>> EEE;
	std::set<std::vector<int>> EET;

	std::set<std::vector<int>> VeVe;
	std::set<std::vector<int>> VeVt;

	std::set<std::vector<int>> EEVe;
	std::set<std::vector<int>> EVeT;
	//std::vector<std::vector<int>> EEET;


	//;

 // all vertices from viable edges
	//std::vector<std::vector<int>> EEVt;
	//
	//std::vector<std::vector<int>> VtTT; // in plane through vertex

	//std::vector<std::vector<int>> P; // plane of primitive (one per facade)
	//std::vector<std::vector<int>> EP;
	//std::vector<std::vector<int>> EEP;

	//std::vector<std::vector<int>> TT; // 

	CombiConfigurations combiS;

	bool findAllEsl() {
		markPotentialSilhouettes();
		markSplitLines();
		storeT();
		storeEE_ET_EVt_VeT_EVe_combis();
	}

	void storeVt() {
		for (Vertex* v : rst->model->vertices) {
			bool wrongside = false;
			for (Primitive* prim : v->triangles) {
				wrongside = false;
				if (glm::dot(rst->window->normal, prim->normal) >= 0) break;
				if (rst->window->intersectsPlaneFromLines(prim->getRayVector())) break;
				wrongside = true;
			}
			if (!wrongside) Vt.insert(v->id);
		}
	}

	void storeT() {
		for (Edge* e : rst->model->edgeVector) {
			bool wrongside = false;
			for (Primitive* prim : e->triangles) {
				wrongside = false;
				if (glm::dot(rst->window->normal, prim->normal) >= 0) break;
				if (rst->window->intersectsPlaneFromLines(prim->getRayVector())) break;
				wrongside = true;
			}
			if (!wrongside) T.insert(e->id);
		}
	}

	void storeEEVe_EVeT_combis() {
		// check for shared vertices??
		for (std::vector<int> vec : EVe) {
			int id1 = vec[0];
			int id2 = vec[1];
			Edge* e1 = rst->model->edgeVector[vec[0]];
			Vertex* v2 = rst->model->vertices[vec[1]];
			for (int id3 : T) {
				Edge* e3 = rst->model->edgeVector[vec[2]];
				if (e3->vertices[0] == v2 || e3->vertices[1] == v2) continue;
				std::vector<glm::dvec3> e1v = { e1->vertices[0]->pos, e1->vertices[1]->pos };
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				Tetrahedron unboundTetra(e1v, e3v);
				if (!unboundTetra.pointInTetra(v2->pos, 1E-7)) continue;

				std::vector<int> E1E3vec;
				if (id1 < id3) E1E3vec = { id1, id3 };
				else E1E3vec = { id3, id1 };

				bool E1T3 = false;
				bool V2T3 = false;
				bool E1E3 = false;
				bool V2E3 = false;

				if (e3->silhouette && id3 > id1) {
					for (std::vector<int> vec2 : EE) {
						if (vec2 == E1E3vec) {
							E1E3 = true;
							break;
						}
					}
					for (std::vector<int> vec2 : EVe) {
						if (vec2 == std::vector<int>{id3, id2}) {
							V2E3 = true;
							break;
						}
					}
					if (E1E3 && V2E3) {
						EEVe.insert({ id1, id2, id3 });
						continue;
					}
				}
				for (std::vector<int> vec2 : ET) {
					if (vec2 == std::vector<int>{id1, id3}) {
						E1T3 = true;
						break;
					}
				}
				for (std::vector<int> vec2 : VeT) {
					if (vec2 == std::vector<int>{id2, id3}) {
						V2T3 = true;
						break;
					}
				}
				if (e1->isInFrontOfEdge(e3, rst->window) && 
					(v2->isInFrontOfVertex(e3->vertices[0], rst->window) || v2->isInFrontOfVertex(e3->vertices[1], rst->window)) && 
					((E1T3 && V2T3) || (E1E3 && V2T3) || (V2E3 && E1T3)))
					EVeT.insert({ id1, id2, id3 });
			}
		}
	}

	void storeEEE_EET_combis() {

		for (std::vector<int> vec : EE) {
			int id1 = vec[0];
			int id2 = vec[1];
			Edge* e1 = rst->model->edgeVector[vec[0]];
			Edge* e2 = rst->model->edgeVector[vec[1]];
			for (int id3 : T) {

				Edge* e3 = rst->model->edgeVector[vec[2]];
				std::vector<glm::dvec3> e1v = { e1->vertices[0]->pos, e1->vertices[1]->pos };
				std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };
				Tetrahedron unboundTetra(e1v, e2v);
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				if (!unboundTetra.segmInTetra(e3v[0], e3v[1], 1E-7)) continue;

				std::vector<int> E1E3vec;
				if (id1 < id3) E1E3vec = { id1, id3 };
				else E1E3vec = { id3, id1 };
				std::vector<int> E2E3vec;
				if (id2 < id3) E2E3vec = { id2, id3 };
				else E2E3vec = { id3, id2 };

				bool E1T3 = false;
				bool E2T3 = false;
				bool E1E3 = false;
				bool E2E3 = false;

				if (e3->silhouette && id3 > id2) {
					for (std::vector<int> vec2 : EE) {
						if (vec2 == E1E3vec) E1E3 = true;
						if (vec2 == E2E3vec) E2E3 = true;
						if (E1E3 && E2E3) break;
					}
					if (E1E3 && E2E3) {
						EEE.insert({ id1, id2, id3 });
						continue;
					}
				}

				for (std::vector<int> vec2 : ET) {
					if (vec2 == std::vector<int>{id1, id3}) E1T3 = true;
					if (vec2 == std::vector<int>{id2, id3}) E2T3 = true;
					if ((E1E3 && E2T3) || (E2E3 && E1T3) || (E1T3 && E2T3)) break;
				}
				if (e1->isInFrontOfEdge(e3, rst->window) && e2->isInFrontOfEdge(e3, rst->window) &&
					((E1T3 && E2T3) || (E1E3 && E2T3) || (E2E3 && E1T3)))
					EET.insert({ id1, id2, id3 });
			}
		}
	}

	void storeEE_ET_EVt_VeT_EVe_combis() {
		for (int T_id1 : T) {
			Edge* e1 = rst->model->edgeVector[T_id1];
			Vertex* e1v1 = e1->vertices[0];
			Vertex* e1v2 = e1->vertices[1];

			for (int T_id2 : T) {
				if (T_id2 <= T_id1) continue;
				Edge* e2 = rst->model->edgeVector[T_id2];
				Vertex* e2v1 = e2->vertices[0];
				Vertex* e2v2 = e2->vertices[1];

				bool e1_for_e2 = false;
				bool e2_for_e1 = false;
				bool e1_for_e2v1 = false;
				bool e1_for_e2v2 = false;
				bool e2_for_e1v1 = false;
				bool e2_for_e1v2 = false;
				bool e2v1_for_e1 = false;
				bool e2v2_for_e1 = false;
				bool e1v1_for_e2 = false;
				bool e1v2_for_e2 = false;

				bool e2v1_for_e1v1 = false;
				bool e2v2_for_e1 = false;
				bool e1v1_for_e2 = false;
				bool e1v2_for_e2 = false;

				if (e1->silhouette) {
					bool side1, side2;
					e1_for_e2v1 = e1->isSilhouetteForPos(e2v1->pos, side1);
					e1_for_e2v2 = e1->isSilhouetteForPos(e2v2->pos, side2);
					e1_for_e2 = e1_for_e2v1 || e1_for_e2v2 || (side1 != side2);
					e1v1_for_e2 = e1_for_e2;
					e1v2_for_e2 = e1_for_e2;
				}
				if (e2->silhouette) {
					bool side1, side2;
					e2_for_e1v1 = e2->isSilhouetteForPos(e1v1->pos, side1);
					e2_for_e1v2 = e2->isSilhouetteForPos(e1v2->pos, side2);
					e2_for_e1 = e2_for_e1v1 || e2_for_e1v2 || (side1 != side2);
					e2v1_for_e1 = e2_for_e1;
					e2v2_for_e1 = e2_for_e1;
				}
				if (!e1_for_e2 && !e2_for_e1 && !e1_for_e2v1 && !e1_for_e2v2 && !e2_for_e1v1 && !e2_for_e1v2) continue;

				Ray e1v1e2v1(e1v1->pos, e2v1->pos);
				Ray e1v1e2v2(e1v1->pos, e2v2->pos);
				Ray e1v2e2v1(e1v2->pos, e2v1->pos);
				Ray e1v2e2v2(e1v2->pos, e2v2->pos);

				if (e1v1 == e2v1) { e1_for_e2v1 = false; e2_for_e1v1 = false; e2v1_for_e1 = false; e1v1_for_e2 = false; }
				if (e1v1 == e2v2) { e1_for_e2v2 = false; e2_for_e1v1 = false; e1v1_for_e2 = false; e2v2_for_e1 = false; }
				if (e1v2 == e2v1) { e2_for_e1v2 = false; e1_for_e2v1 = false; e1v2_for_e2 = false; e2v1_for_e1 = false; }
				if (e1v2 == e2v2) { e2_for_e1v2 = false; e1_for_e2v2 = false; e1v2_for_e2 = false; e2v2_for_e1 = false; }

				if (!e1_for_e2 && !e2_for_e1 && !e1_for_e2v1 && !e1_for_e2v2 && !e2_for_e1v1 && !e2_for_e1v2) continue;

				if (e1_for_e2 || e2_for_e1) {
					std::vector<Ray> lines = { e1v1e2v1, e1v1e2v2, e1v2e2v1, e1v2e2v2};
					if (!rst->window->intersectsPlaneFromLines(lines)) continue;
				}
				if (e1_for_e2v1 || e2v1_for_e1) {
					std::vector<Ray> lines = { e1v1e2v1, e1v2e2v1 };
					if (!rst->window->intersectsPlaneFromLines(lines)) {
						e1_for_e2v1 = false;
						e2v1_for_e1 = false;
					}
				}
				if (e1_for_e2v2 || e2v2_for_e1) {
					std::vector<Ray> lines = { e1v1e2v2, e1v2e2v2 };
					if (!rst->window->intersectsPlaneFromLines(lines)) {
						e1_for_e2v1 = false;
						e2v2_for_e1 = false;
					}
				}
				if (e2_for_e1v1 || e1v1_for_e2) {
					std::vector<Ray> lines = { e1v1e2v1, e1v1e2v2 };
					if (!rst->window->intersectsPlaneFromLines(lines)) {
						e1_for_e2v1 = false;
						e1v1_for_e2 = false;
					}
				}
				if (e2_for_e1v2 || e1v2_for_e2) {
					std::vector<Ray> lines = { e1v2e2v1, e1v2e2v2 };
					if (!rst->window->intersectsPlaneFromLines(lines)) {
						e2_for_e1v2 = false;
						e1v2_for_e2 = false;
					}
				}
				if (!e1_for_e2 && !e2_for_e1 && !e1_for_e2v1 && !e1_for_e2v2 && !e2_for_e1v1 && !e2_for_e1v2) continue;

				if (e1_for_e2 && e2_for_e1) {
					EE.insert({ T_id1, T_id2 });
					if (e1_for_e2v1) EVe.insert({ e1->id, e2v1->id });
					if (e1_for_e2v2) EVe.insert({ e1->id, e2v2->id });
					if (e2_for_e1v1) EVe.insert({ e2->id, e1v1->id });
					if (e2_for_e1v2) EVe.insert({ e2->id, e1v2->id });
					if (rst->window->inBounds(e1v2e2v2)) {
						if (e2_for_e1v2 && e1_for_e2v2) VeVe.insert({ e1v2->id, e2v2->id });
						else if (e2_for_e1v2 && e2v2->isInFrontOfVertex(e1v2, rst->window)) VeVt.insert({ e2v2->id, e1v2->id });
						else if (e1_for_e2v2 && e1v2->isInFrontOfVertex(e2v2, rst->window)) VeVt.insert({ e1v2->id, e2v2->id });
					}
					if (rst->window->inBounds(e1v1e2v2)) {
						if (e2_for_e1v1 && e1_for_e2v2) VeVe.insert({ e1v1->id, e2v2->id });
						else if (e2_for_e1v1 && e2v2->isInFrontOfVertex(e1v1, rst->window)) VeVt.insert({ e2v2->id, e1v1->id });
						else if (e1_for_e2v2 && e1v1->isInFrontOfVertex(e2v2, rst->window)) VeVt.insert({ e1v1->id, e2v2->id });
					}
					if (rst->window->inBounds(e1v2e2v1)) {
						if (e2_for_e1v2 && e1_for_e2v1) VeVe.insert({ e1v2->id, e2v1->id });
						else if (e2_for_e1v2 && e2v1->isInFrontOfVertex(e1v2, rst->window)) VeVt.insert({ e2v1->id, e1v2->id });
						else if (e1_for_e2v1 && e1v2->isInFrontOfVertex(e2v1, rst->window)) VeVt.insert({ e1v2->id, e2v1->id });
					}
					if (rst->window->inBounds(e1v1e2v1)) {
						if (e2_for_e1v1 && e1_for_e2v1) VeVe.insert({ e1v1->id, e2v1->id });
						else if (e2_for_e1v1 && e2v1->isInFrontOfVertex(e1v1, rst->window)) VeVt.insert({ e2v1->id, e1v1->id });
						else if (e1_for_e2v1 && e1v1->isInFrontOfVertex(e2v1, rst->window)) VeVt.insert({ e1v1->id, e2v1->id });
					}
				}
				else if (e1_for_e2 && e1->isInFrontOfEdge(e2, rst->window)) {
					ET.insert({ T_id1, T_id2 });
					if (e1_for_e2v1 && e1->isInFrontOfVertex(e2v1, rst->window)) {
						EVt.insert({ e1->id, e2v1->id });
						if (e1v1_for_e2 && rst->window->inBounds(e1v1e2v1) && e1v1->isInFrontOfVertex(e2v1, rst->window)) 
							VeVt.insert({ e1v1->id, e2v1->id });
						if (e1v2_for_e2 && rst->window->inBounds(e1v2e2v1) && e1v2->isInFrontOfVertex(e2v1, rst->window))
							VeVt.insert({ e1v2->id, e2v1->id });
					}
					if (e1_for_e2v2 && e1->isInFrontOfVertex(e2v2, rst->window)) {
						EVt.insert({ e1->id, e2v2->id });
						if (e1v1_for_e2 && rst->window->inBounds(e1v1e2v2) && e1v1->isInFrontOfVertex(e2v2, rst->window))
							VeVt.insert({ e1v1->id, e2v2->id });
						if (e1v2_for_e2 && rst->window->inBounds(e1v2e2v2) && e1v2->isInFrontOfVertex(e2v2, rst->window))
							VeVt.insert({ e1v2->id, e2v2->id });
					}
					if (e1v1_for_e2 && (e1v1->isInFrontOfVertex(e2v1, rst->window) || e1v1->isInFrontOfVertex(e2v2, rst->window)))
						VeT.insert({ e1v1->id, e2->id });
					if (e1v2_for_e2 && (e1v2->isInFrontOfVertex(e2v1, rst->window) || e1v2->isInFrontOfVertex(e2v2, rst->window)))
						VeT.insert({ e1v2->id, e2->id });
				}
				else if (e2_for_e1 && e2->isInFrontOfEdge(e1, rst->window)) {
					ET.insert({ T_id2, T_id1 });
					if (e2_for_e1v1 && e2->isInFrontOfVertex(e1v1, rst->window)) {
						EVt.insert({ e2->id, e1v1->id });
						if (e2v1_for_e1 && rst->window->inBounds(e1v1e2v1) && e2v1->isInFrontOfVertex(e1v1, rst->window))
							VeVt.insert({ e2v1->id, e1v1->id });
						if (e2v2_for_e1 && rst->window->inBounds(e1v1e2v2) && e2v2->isInFrontOfVertex(e1v1, rst->window))
							VeVt.insert({ e2v2->id, e1v1->id });
					}
					if (e2_for_e1v2 && e2->isInFrontOfVertex(e1v2, rst->window)) {
						EVt.insert({ e2->id, e1v2->id });
						if (e2v1_for_e1 && rst->window->inBounds(e1v2e2v1) && e2v1->isInFrontOfVertex(e1v2, rst->window))
							VeVt.insert({ e2v1->id, e1v2->id });
						if (e2v2_for_e1 && rst->window->inBounds(e1v2e2v2) && e2v2->isInFrontOfVertex(e1v2, rst->window))
							VeVt.insert({ e2v2->id, e1v2->id });
					}
					if (e2v1_for_e1 && (e2v1->isInFrontOfVertex(e1v1, rst->window) || e2v1->isInFrontOfVertex(e1v2, rst->window)))
						VeT.insert({ e2v1->id, e1->id });
					if (e2v2_for_e1 && (e2v2->isInFrontOfVertex(e1v1, rst->window) || e2v2->isInFrontOfVertex(e1v2, rst->window)))
						VeT.insert({ e2v2->id, e1->id });
				}
			}
		}
		// remove the EVt instances that are also in EVe as they are misclassified
		std::vector<std::vector<int>> remove;
		for (std::vector<int> EVtvec : EVt) {
			for (std::vector<int> EVevec : EVe) {
				if (EVtvec == EVevec) remove.push_back(EVtvec);
			}
		}
		for (std::vector<int> rm : remove) EVt.erase(rm);
		// same for VeVt also in VeVe
		std::vector<std::vector<int>> remove2;
		for (std::vector<int> VeVtvec : VeVt) {
			for (std::vector<int> VeVevec : VeVe) {
				if (VeVtvec == VeVevec) remove.push_back(VeVevec);
			}
		}
		for (std::vector<int> rm : remove2) VeVt.erase(rm);
	}

	bool markSplitLines() {
		for (Split split : rst->splitters) {
			split.edge->splitline = true;
			for (Vertex* v : split.edge->vertices) v->splitline = true;
		}
	}

	bool markPotentialSilhouettes() {
		for (Edge* e : rst->model->edgeVector) {
			if (!e->convex()) continue;
			int count = 0;
			for (glm::dvec3& point : rst->window->vertices) {
				bool side;
				if (e->isSilhouetteForPos(point, side)) {
					e->silhouette = true;
					break;
				}
				if (side) count++;
			}
			if (count > 0 && count <rst->window->vertices.size()) e->silhouette = true;
			if (e->silhouette) {
				for (Vertex* v : e->vertices) v->silhouette = true;
			}
		}
	}

	bool checkESL(ESLCandidate& esl) {

	}

	// should store potentially two lines?? for all ESLs
	bool checkRaysThrough4Lines(ESLCandidate& esl, bool vischeck) {

		std::vector<Line4> intersectLines = Lines4Finder::find(esl.lines4);
		for (int i = 0; i < intersectLines.size(); i++) {
			ESLCandidate esltry = esl;
			esltry.inPlane = false;
			intersectLines[i].get3DfromPlucker();
			esltry.putRay(intersectLines[i]);
			if (checkESL(esltry)) return true;
		}
		return false;
	}



};