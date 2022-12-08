#pragma once
#include <chrono>
#include "cache.h"
#include "model.h"
#include "node.h"
#include "findLineThroughFour.h"
#include "combinations.h"
#include "rst.h"
#include "eslChecker.h"
#include "eslCandidate.h"
#include "extremalStabbingLine.h"
#include "splitSide.h"
#include "combiConfigurations.h"
#include <unordered_set>
#include "eslCandidatePre.h"

class ESLFindAll {
public:
	//std::vector<ESLCandidate> eslsNoVis;
	//std::vector<ESLCandidate> esls;
	//std::vector<Edge*> silhouetteEdges;
	bool checkHardCombis = false;
	std::vector<std::vector<ExtremalStabbingLine>> eslsPerPrimitive;
	std::vector<std::vector<Ray>> raysPerPrimitive;
	std::vector<ESLCandidatePre> esls;

	ESLFindAll() {};

	ESLFindAll(RaySpaceTree* rst) :rst(rst)
	{
		//eslChecker = ESLChecker(node, prim, rst, printAll);
	};

	void find() {
		findAllEsl();
	};

private:
	
	struct hash {
		std::size_t operator()(std::vector<uint32_t> const& vec) const {
			std::size_t seed = vec.size();
			for (auto x : vec) {
				x = ((x >> 16) ^ x) * 0x45d9f3b;
				x = ((x >> 16) ^ x) * 0x45d9f3b;
				x = (x >> 16) ^ x;
				seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};


	RaySpaceTree* rst;

	ESLChecker eslChecker;

	std::vector<Vertex*> silhouetteVertices;
	std::vector<SplitSide> splitLines;

	std::set<int> T; // all viable edges (can be seen from window, depends on planes of primitives attached)
	std::set<int> Vt;
	std::unordered_set<std::vector<uint32_t>, hash> ET;
	std::unordered_set<std::vector<uint32_t>, hash> EE;
	std::unordered_set<std::vector<uint32_t>, hash> VeT;
	std::unordered_set<std::vector<uint32_t>, hash> EVe;
	std::unordered_set<std::vector<uint32_t>, hash> EVt;
	std::unordered_set<std::vector<uint32_t>, hash> EEE;
	std::unordered_set<std::vector<uint32_t>, hash> VeVe;
	std::unordered_set<std::vector<uint32_t>, hash> VeVt;
	std::unordered_set<std::vector<uint32_t>, hash> EEVe;
	std::unordered_set<std::vector<uint32_t>, hash> EEEE;

	std::unordered_set<std::vector<uint32_t>, hash> WE;
	std::unordered_set<std::vector<uint32_t>, hash> WEE;
	std::unordered_set<std::vector<uint32_t>, hash> WEEE;
	std::unordered_set<std::vector<uint32_t>, hash> WVe;
	std::unordered_set<std::vector<uint32_t>, hash> WEVe;
	std::unordered_set<std::vector<uint32_t>, hash> VwVe;
	std::unordered_set<std::vector<uint32_t>, hash> VwE;
	std::unordered_set<std::vector<uint32_t>, hash> VwEE;

	// all vertices from viable edges
	//std::vector<std::vector<uint32_t>> VtTT; // in plane through vertex
	//std::vector<std::vector<uint32_t>> P; // plane of primitive (one per facade)
	//std::vector<std::vector<uint32_t>> EP;
	//std::vector<std::vector<uint32_t>> EEP;
	//std::vector<std::vector<uint32_t>> TT; // 

	//std::unordered_set<std::vector<uint32_t>, hash> EET;
	//std::unordered_set<std::vector<uint32_t>, hash> EVeT;
	//std::unordered_set<std::vector<uint32_t>, hash> EEVt;
	//std::unordered_set<std::vector<uint32_t>, hash> EEET;
	//std::unordered_set<std::vector<uint32_t>, hash> WT;
	//std::unordered_set<std::vector<uint32_t>, hash> WET;
	//std::unordered_set<std::vector<uint32_t>, hash> WEET;
	//std::unordered_set<std::vector<uint32_t>, hash> WVt;
	//std::unordered_set<std::vector<uint32_t>, hash> WEVt;
	//std::unordered_set<std::vector<uint32_t>, hash> WVeT;
	//std::unordered_set<std::vector<uint32_t>, hash> VwVt;
	//std::unordered_set<std::vector<uint32_t>, hash> VwT;
	//std::unordered_set<std::vector<uint32_t>, hash> VwET;

	CombiConfigurations combiS;

	bool findAllEsl() {

		auto start_time = std::chrono::high_resolution_clock::now();

		// remove veve if same primitive - but add another check that does this? (also add window check?)
		// also for edge/edge edge/vertex combis?
		std::cout << "Nr of Primitives: " << rst->model->primsize << std::endl;

		markPotentialSilhouettes();
		markSplitLines();

		storeT();											// SSST SSSE
		storeVt();											// SSVt SSVe
		std::cout << std::endl << "EE_ET_EVe_EVt_VeT_VeVt_VeVe " << T.size() << "x" << T.size() << std::endl;
		storeEE_ET_EVe_EVt_VeT_VeVt_VeVe_combis();			// SSEE SSET SEVt SVeT SEVe VeVe VeVt
		std::cout << std::endl << "EEE_EET " << EE.size() << "x" << T.size() << std::endl;
		storeEEE_EET_combis();								// SEEE SEET
		std::cout << std::endl << "EEVe_EVeT " << EVe.size() << "x" << T.size() << std::endl;
		storeEEVe_EVeT_combis();							// EEVe EVeT
		//std::cout << std::endl << "EEVt" << EE.size() << "x" << Vt.size() << std::endl;
		//storeEEVt_combis();									// EEVt
		std::cout << std::endl << "EEEE_EEET " << EEE.size() << "x" << T.size() << std::endl;
		storeEEEE_EEET_combis();							// EEEE EEET
		std::cout << std::endl << "WE_WT_WVt_WVe " << T.size() << "x" << 4 << std::endl;
		storeWE_WT_WVt_WVe_combis();						// WSSE WSST WSVe WSVt
		std::cout << std::endl << "VwVe_VwVt_VwE_VwT " << T.size() << "x" << 4 << std::endl;
		storeVwVe_VwVt_VwE_VwTcombis();						// VwVe VwVt VwST VwSE
		std::cout << std::endl << "VwEE_VwET " << VwE.size() << "x" << T.size() << std::endl;
		storeVwEE_VwETcombis();								// VwEE VwET
		std::cout << std::endl << "WEVt_WEVe " << WE.size() << "x" << Vt.size() << std::endl;
		storeWEVt_WEVe_combis();							// WEVt WEve
		//std::cout << std::endl << "WVeT " << WVe.size() << "x" << T.size() << std::endl;
		//storeWVeTcombis();									// WVeT
		std::cout << std::endl << "WEE_WET " << WE.size() << "x" << T.size() << std::endl;
		storeWEE_WET_combis();								// WSEE	WSET
		std::cout << std::endl << "WEEE_WEET " << WEE.size() << "x" << T.size() << std::endl;
		storeWEEE_WEET_combis();							// WEEE WEET		

		// check if inplane stuff is covered! or add seperately
		// or make sure that if edgeA is in plane of 1 of 2 silh triangles of edgeB, edgeB is silh for edgeA
		// especially extra check to see if part of same primitive -> definitely in plane
		makeESLsFromCombi();

		auto end_time = std::chrono::high_resolution_clock::now();

		std::cout << std::endl << "Found " << esls.size() << " ESLs in " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms!" << std::endl;
		//for (ESLCandidatePre &esl : esls) {
		//	std::cout << esl.Ve[0]->id << " " << esl.Ve[1]->id << std::endl;
		//}

		return true;
	}

	void fillPrimsFromESLs(ESLCandidatePre& candidate) {
		candidate.findPrimsForESLs();
		bool addESL = false;
		for (Line4& r : candidate.rays) {
			if (r.prims.size()) addESL = true;
			for (int prim : r.prims) raysPerPrimitive[prim].push_back(r);
		}
		if (addESL) {
			candidate.id = esls.size();
			esls.push_back(candidate);
		}
	}

	void makeESLsFromCombi() {
		raysPerPrimitive.resize(rst->model->primsize);
		for (std::vector<uint32_t> vwve : VwVe) fillPrimsFromESLs(ESLCandidatePre(vwve, 2, rst, 1, 0, 1, 0));
		for (std::vector<uint32_t> veve : VeVe) fillPrimsFromESLs(ESLCandidatePre(veve, 2, rst, 0, 0, 2, 0));
		for (std::vector<uint32_t> vwee : VwEE) fillPrimsFromESLs(ESLCandidatePre(vwee, 1, rst, 1, 0, 0, 2));
		for (std::vector<uint32_t> weee : WEEE) fillPrimsFromESLs(ESLCandidatePre(weee, 0, rst, 0, 1, 0, 3));
		for (std::vector<uint32_t> weve : WEVe) fillPrimsFromESLs(ESLCandidatePre(weve, 1, rst, 0, 1, 1, 1));
		for (std::vector<uint32_t> eeve : EEVe) fillPrimsFromESLs(ESLCandidatePre(eeve, 1, rst, 0, 0, 1, 2));
		for (std::vector<uint32_t> eeee : EEEE) fillPrimsFromESLs(ESLCandidatePre(eeee, 0, rst, 0, 0, 0, 4));
	}

	bool markSplitLines() {
		for (Split split : rst->splitters) {
			split.edge->splitline = true;
			for (Vertex* v : split.edge->vertices) v->splitline = true;
		}
		return true;
	}

	bool potentialSilhouetteForWindow(Edge* e) {
		if (!e->convex()) return false;

		int count = 0;
		for (glm::dvec3& point : rst->window->vertices) {
			bool side;
			if (e->isSilhouetteForPos(point, side)) return true;
			if (side) count++;;
		}
		if (count == 0 || count == rst->window->vertices.size()) return false;
		return true;
	}

	bool markPotentialSilhouettes() {
		for (Edge* e : rst->model->edgeVector) {
			e->silhouette = potentialSilhouetteForWindow(e);
		}
		for (Vertex* v : rst->model->vertices) {
			int count = 0;
			for (Edge* e : v->edges) {
				if (e->silhouette) count++;
			}
			if (count >= 2) v->silhouette = true;
		}
		return true;
	}

	void storeVt() {
		for (Vertex* v : rst->model->vertices) {
			for (Primitive* prim : v->triangles) {
				if (glm::dot(rst->window->normal, prim->normal) > 0 ||
					rst->window->planeIntersectionFromLines(prim->getRayVector())) {
					Vt.insert(v->id);
					break;
				}
			}
		}
	}

	void storeT() {
		for (Edge* e : rst->model->edgeVector) {
			for (Primitive* prim : e->triangles) {
				if (glm::dot(rst->window->normal, prim->normal) > 0 ||
					rst->window->planeIntersectionFromLines(prim->getRayVector())) {
					T.insert(e->id);
					break;
				}
			}
		}
	}

	//void storeWVeTcombis() {
	//	int count = 0;
	//	for (std::vector<uint32_t> vec : WVe) {
	//		if (count % 1000 == 0) std::cout << count << " ";
	//
	//		uint32_t id1 = vec[0];
	//		uint32_t id2 = vec[1];
	//		std::vector<glm::dvec3> w1v = { rst->window->vertices[id1], rst->window->vertices[(id1 + 1) % rst->window->vertices.size()] };
	//		Vertex* v2 = rst->model->vertices[id2];
	//		TriWedge w1v2(v2->pos, w1v[0], w1v[1], glm::vec3(0, 0, 0));
	//
	//		for (uint32_t id3 : T) {
	//			Edge* e3 = rst->model->edgeVector[id3];
	//			std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
	//			glm::dvec3 pt;
	//			if (!w1v2.raySegmentIntersection(e3v[0], e3v[1], pt)) continue;
	//			if (WEVe.count( { id1, id3, id2 })) continue;
	//			bool W1T3 = WT.count({id1, id3});
	//			bool Ve2T3 = VeT.count({id2, id3});
	//			bool W1E3 = WE.count({ id1, id3 });
	//			bool E3Ve2 = EVe.count({ id3, id2 });
	//			if (rst->window->isInFrontOf(v2->pos, e3v) && (W1T3 || W1E3) && (Ve2T3 || E3Ve2)) WVeT.insert({ id1, id2, id3 });
	//		}
	//	}
	//}

	void storeWEVt_WEVe_combis() {
		int count = 0;
		for (std::vector<uint32_t> vec : WE) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e2 = rst->model->edgeVector[id2];

			std::vector<glm::dvec3> w1v = rst->window->vertexEdges[id1];
			std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };

			for (uint32_t id3 : Vt) {
				Vertex* v3 = rst->model->vertices[id3];

				TriWedge e1v3(v3->pos, w1v[0], w1v[1], glm::vec3(0, 0, 0));
				glm::dvec3 pt;
				if (!e1v3.raySegmentIntersection(e2v[0], e2v[1], pt)) continue;

				bool W1Ve3 = false;
				bool E2Ve3 = false;
				if (v3->silhouette) {
					W1Ve3 = WVe.count({id1, id3});
					E2Ve3 = EVe.count({id2, id3});
					if (W1Ve3 && E2Ve3) {
						WEVe.insert({ id1, id2, id3 });
						continue;
					}
				}
				//if (!rst->window->isInFrontOf(e2v, v3->pos)) continue;
				//bool W1Vt3 = WVt.count({id1, id3});
				//bool E2Vt3 = EVt.count({id2, id3});
				//if ((W1Vt3 || W1Ve3) && (E2Ve3 || E2Vt3)) WEVt.insert({ id1, id2, id3 });
			}
		}
	}

	void storeVwEE_VwETcombis() {
		int count = 0;
		for (std::vector<uint32_t> vec : VwE) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			glm::dvec3 v1 = rst->window->vertices[id1];
			Edge* e2 = rst->model->edgeVector[id2];
			std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };
			TriWedge e1v2(v1, e2v[0], e2v[1], glm::vec3(0, 0, 0));
			for (uint32_t id3 : T) {
				Edge* e3 = rst->model->edgeVector[id3];
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				glm::dvec3 pt;
				if (!e1v2.raySegmentIntersection(e3v[0], e3v[1], pt)) continue;

				bool E2E3 = false;
				bool Vw1E3 = false;
				if (e3->silhouette && id3 > id2) {
					E2E3 = EE.count({ id2, id3 });
					Vw1E3 = VwE.count({ id1, id3 });
					if (E2E3 && Vw1E3) {
						VwEE.insert({id1, id2, id3});
						continue;
					}
				}
				//if (!rst->window->isInFrontOf(e2v, e3v)) continue;
				//bool E2T3 = false;
				//bool Vw1T3 = false;
				//E2T3 = ET.count({ id2, id3 });
				//Vw1T3 = VwT.count({ id1, id3 });
				//if ((E2T3 || E2E3) && (Vw1E3 || Vw1T3)) VwET.insert({ id1, id2, id3 });
			}
		}
	}
	void storeWEE_WET_combis() {
		int count = 0;
		for (std::vector<uint32_t> vec : WE) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e2 = rst->model->edgeVector[id2];

			std::vector<glm::dvec3> w1v = rst->window->vertexEdges[id1];
			std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };
			Tetrahedron unboundTetra(w1v, e2v);

			for (uint32_t id3 : T) {
				Edge* e3 = rst->model->edgeVector[id3];
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				if (!unboundTetra.segmInTetra(e3v[0], e3v[1], 1E-7)) continue;

				//bool W1T3 = false;
				//bool E2T3 = false;
				bool W1E3 = false;
				bool E2E3 = false;

				if (e3->silhouette && id3 > id2) {
					E2E3 = EE.count({id2, id3});
					W1E3 = WE.count({id1, id3});
					if (W1E3 && E2E3) {
						WEE.insert({ id1, id2, id3 });
						continue;
					}
				}
				//if (!rst->window->isInFrontOf(e2v, e3v)) continue;
				//W1T3 = WT.count({id1, id3});
				//E2T3 = ET.count({id2, id3});
				//if ((W1T3 && E2T3) || (W1T3 && E2E3) || (W1E3 && E2T3)) WET.insert({ id1, id2, id3 });
			}
		}
	}

	void storeWEEE_WEET_combis() {
		int count = 0;
		for (std::vector<uint32_t> vec : WEE) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			uint32_t id3 = vec[2];
			Edge* e2 = rst->model->edgeVector[id2];
			Edge* e3 = rst->model->edgeVector[id3];

			for (uint32_t id4 : T) {
				Edge* e4 = rst->model->edgeVector[id4];

				bool W1E2E4 = false;
				bool W1E3E4 = false;
				bool E2E3E4 = false;

				if (e4->silhouette && id4 > id3) {
					W1E2E4 = WEE.count({id1, id2, id4});
					W1E3E4 = WEE.count({id1, id3, id4});
					E2E3E4 = EEE.count({id2, id3, id4});
					if (W1E2E4 && W1E3E4 && E2E3E4) {
						WEEE.insert({ id1, id2, id3, id4 });
						continue;
					}
				}
				//if (!rst->window->isInFrontOf(e2->getVecPos(), e4->getVecPos()) || !rst->window->isInFrontOf(e3->getVecPos(), e4->getVecPos())) continue;
				//bool W1E2T4 = WT.count({id1, id2, id4});
				//bool W1E3T4 = ET.count({id1, id3, id4});
				//bool E2E3T4 = EET.count({id2, id3, id4});
				//if ((E2E3T4 || E2E3E4) && (W1E2E4 || W1E2T4) && (W1E3E4 || W1E3T4)) WEET.insert({ id1, id2, id3 });
			}
		}
	}

	void storeVwVe_VwVt_VwE_VwTcombis() {
		int count = 0;
		for (uint32_t id1 : Vt) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			Vertex* v2 = rst->model->vertices[id1];
			for (uint32_t i = 0; i < rst->window->vertices.size(); i++) {
				bool W1aVe2 = false;
				bool W1bVe2 = false;
				uint32_t j = (i + 1) % rst->window->vertices.size();
				if (v2->silhouette) {
					W1aVe2 = WVe.count({ i, id1 });
					W1bVe2 = WVe.count({ j, id1 });
					if (W1aVe2 && W1bVe2) {
						VwVe.insert({ i, id1 });
						continue;
					}
				}
				//bool W1aVt2 = false;
				//bool W1bVt2 = false;
				//W1aVt2 = WVt.count({ i, id1 });
				//W1bVt2 = WVt.count({ j, id1 });
				//if ((W1aVe2 || W1aVt2) && (W1bVe2 || W1bVt2)) VwVt.insert({ i, id1 });

			}
		}
		for (uint32_t id1 : T) {
			Edge* e = rst->model->edgeVector[id1];
			for (uint32_t i = 0; i < rst->window->vertices.size(); i++) {
				bool W1aE2 = false;
				bool W1bE2 = false;
				uint32_t j = (i + 1) % rst->window->vertices.size();
				if (e->silhouette) {
					W1aE2 = WE.count({ i, id1 });
					W1bE2 = WE.count({ j, id1 });
					if (W1aE2 && W1bE2) {
						VwE.insert({ i, id1 });
						continue;
					}
				}
				//bool W1aT2 = false;
				//bool W1bT2 = false;
				//W1aT2 = WT.count({ i, id1 });
				//W1bT2 = WT.count({ j, id1 });
				//if ((W1aE2 || W1aT2) && (W1bE2 || W1bT2)) VwT.insert({ i, id1 });
			}
		}
	}

	void storeWE_WT_WVt_WVe_combis() {
		int count = 0;
		for (uint32_t id1 : T) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			Edge* e1 = rst->model->edgeVector[id1];

			for (uint32_t i = 0; i < rst->window->vertices.size(); i++) {
				std::vector<glm::dvec3> w1 = rst->window->vertexEdges[i];

				if (e1->silhouette && e1->isSilhouetteForPrim(w1)) {
					WE.insert({ i, id1 });
					WVe.insert({i, uint32_t(e1->vertices[0]->id)});
					WVe.insert({ i, uint32_t(e1->vertices[1]->id) });
				}
				//else {
				//	WT.insert({ i, id1 });
				//	WVt.insert({ i, uint32_t(e1->vertices[0]->id) });
				//	WVt.insert({ i, uint32_t(e1->vertices[1]->id) });
				//}
			}
		}
		//std::vector<std::vector<uint32_t>> remove;
		//for (std::vector<uint32_t> WVtvec : WVt) {
		//	if (WVe.count(WVtvec)) remove.push_back(WVtvec);
		//}
		//for (std::vector<uint32_t> rem : remove) WVt.erase(rem);
	}

	void storeEEEE_EEET_combis() {
		int count = 0;
		for (std::vector<uint32_t> vec : EEE) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			uint32_t id3 = vec[2];
			Edge* e1 = rst->model->edgeVector[id1];
			Edge* e2 = rst->model->edgeVector[id2];
			Edge* e3 = rst->model->edgeVector[id3];

			for (uint32_t id4 : T) {
				if (id1 == id4 || id2 == id4 || id3 == id4) continue;
				Edge* e4 = rst->model->edgeVector[id4];

				bool E1E2E4 = false;
				bool E2E3E4 = false;
				bool E1E3E4 = false;

				if (e4->silhouette && id4 > id3) {
					E1E2E4 = EEE.count({ id1, id2, id4 });
					E2E3E4 = EEE.count({ id2, id3, id4 });
					E1E3E4 = EEE.count({ id1, id3, id4 });
					if (E1E2E4 && E2E3E4 && E1E3E4) {
						EEEE.insert({ id1, id2, id3, id4 });
						continue;
					}
				}

				//if (!(rst->window->isInFrontOf(e1->getVecPos(), e4->getVecPos()) && rst->window->isInFrontOf(e2->getVecPos(), e4->getVecPos()) && rst->window->isInFrontOf(e3->getVecPos(), e4->getVecPos()))) continue;

				//bool E1E2T4 = EET.count({ id1, id2, id4 });
				//if (!E1E2T4 && !E1E2E4) continue;
				//bool E2E3T4 = EET.count({ id2, id3, id4 });
				//if (!E2E3T4 && !E2E3E4) continue;
				//bool E1E3T4 = EET.count({ id1, id3, id4 });
				//if (!E1E3E4 && !E1E3T4) continue;
				//EEET.insert({ id1, id2, id3, id4 });
			}
		}
	}

	//void storeEEVt_combis() {
	//	int count = 0;
	//	for (std::vector<uint32_t> vec : EE) {
	//		if (count % 1000 == 0) std::cout << count << " "; count++;
	//
	//		uint32_t id1 = vec[0];
	//		uint32_t id2 = vec[1];
	//		Edge* e1 = rst->model->edgeVector[vec[0]];
	//		Edge* e2 = rst->model->edgeVector[vec[1]];
	//		std::vector<glm::dvec3> e1v = e1->getVecPos();
	//		std::vector<glm::dvec3> e2v = e2->getVecPos();
	//
	//		for (uint32_t id3 : Vt) {
	//			Vertex* v3 = rst->model->vertices[id3];
	//
	//			if (e1->vertices[0] == v3 || e1->vertices[1] == v3 ||
	//				e2->vertices[0] == v3 || e2->vertices[1]) continue;
	//
	//			TriWedge e1v3(v3->pos, e1v[0], e1v[1], glm::vec3(0, 0, 0));
	//			glm::dvec3 pt;
	//			if (!e1v3.raySegmentIntersection(e2v[0], e2v[1], pt)) continue;
	//
	//			if (EEVe.count({ id1, id2, id3 })) continue;
	//
	//			bool E1Vt3 = EVt.count({ id1, id3 });
	//			bool E2Vt3 = EVt.count({ id2, id3 });
	//			if (E1Vt3 && E2Vt3) {
	//				EEVt.insert({ id1, id2, id3 });
	//				continue;
	//			}
	//			bool E1Ve3 = EVe.count({ id1, id3 });
	//			bool E2Ve3 = EVe.count({ id2, id3 });
	//			if (!rst->window->isInFrontOf(e1v, v3->pos) || !rst->window->isInFrontOf(e2v, v3->pos)) continue;
	//			if ((E1Vt3 || E1Ve3) && (E2Vt3 || E2Ve3)) EEVt.insert({ id1, id2, id3 });
	//		}
	//	}
	//}

	void storeEEVe_EVeT_combis() {
		// check for shared vertices??
		int count = 0;
		for (std::vector<uint32_t> vec : EVe) {
			if (count % 1000 == 0) std::cout << count << " "; count++;

			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e1 = rst->model->edgeVector[id1];
			Vertex* v2 = rst->model->vertices[id2];
			std::vector<glm::dvec3> e1v = { e1->vertices[0]->pos, e1->vertices[1]->pos };
			TriWedge e1v2(v2->pos, e1v[0], e1v[1], glm::vec3(0, 0, 0));

			for (uint32_t id3 : T) {
				Edge* e3 = rst->model->edgeVector[id3];
				if (e3->vertices[0] == v2 || e3->vertices[1] == v2) continue;
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				glm::dvec3 pt;
				if (!e1v2.raySegmentIntersection(e3v[0], e3v[1], pt)) continue;

				//bool E1T3 = false;
				//bool V2T3 = false;
				bool E1E3 = false;
				bool V2E3 = false;

				if (e3->silhouette && id3 > id1) {
					E1E3 = EE.count({ id1, id3 });
					V2E3 = EVe.count({ id3, id2 });
					if (E1E3 && V2E3) {
						EEVe.insert({ id1, id3, id2 });
						continue;
					}
				}

				//if (rst->window->isInFrontOf(e1v, e3v) && rst->window->isInFrontOf(v2->pos, e3v)) {
				//	E1T3 = ET.count({ id1, id3 });
				//	V2T3 = VeT.count({ id2, id3 });
				//	if ((E1T3 || E1E3) && (V2E3 || V2T3)) EVeT.insert({ id1, id2, id3 });					
				//}
			}
		}
	}

	void storeEEE_EET_combis() {
		int count = 0;
		for (std::vector<uint32_t> vec : EE) {
			if (count % 1000 == 0) std::cout << count << " "; count++;
			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e1 = rst->model->edgeVector[id1];
			Edge* e2 = rst->model->edgeVector[id2];
			std::vector<glm::dvec3> e1v = { e1->vertices[0]->pos, e1->vertices[1]->pos };
			std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };
			Tetrahedron unboundTetra(e1v, e2v);

			for (uint32_t id3 : T) {
				if (id1 == id3 || id2 == id3) continue;

				Edge* e3 = rst->model->edgeVector[id3];
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				if (!unboundTetra.segmInTetra(e3v[0], e3v[1], 1E-7)) continue;


				bool E1E3 = false;
				bool E2E3 = false;

				if (e3->silhouette && id3 > id2) {
					E1E3 = EE.count({ id1, id3 });
					E2E3 = EE.count({ id2, id3 });
					if (E1E3 && E2E3) {
						EEE.insert({ id1, id2, id3 });
						continue;
					}
				}

				//if (!(rst->window->isInFrontOf(e1v, e3v) && rst->window->isInFrontOf(e2v, e3v))) continue;
				//bool E1T3 = ET.count({ id1, id3 });
				//bool E2T3 = ET.count({ id2, id3 });
				//if ((E1T3 || E1E3) && (E2T3 || E2E3)) EET.insert({ id1, id2, id3 });
			}
		}
	}

	bool vectorAndCheck(std::vector<bool>& vec) {
		for (bool x : vec) if (!x) return false;
		return true;
	}

	bool vectorOrCheck(std::vector<bool>& vec) {
		for (bool x : vec) if (x) return true;
		return false;
	}

	bool vectorAndCheck(std::vector<std::vector<bool>>& vec) {
		for (std::vector<bool> x : vec) if (!vectorAndCheck(x)) return false;
		return true;
	}

	bool vectorOrCheck(std::vector<std::vector<bool>>& vec) {
		for (std::vector<bool> x : vec) if (vectorAndCheck(x)) return true;
		return false;
	}

	template <class T>
	std::vector<T> flatten(std::vector<std::vector<T>>& vec) {
		std::vector<T> flat;
		for (std::vector<T> v : vec) {
			for (T t : v) {
				flat.push_back(t);
			}
		}
		return flat;
	}

	template <class T>
	std::vector<T> getCol(std::vector<std::vector<T>>& vec, int i) {
		std::vector<T> slice;
		for (std::vector<T> v : vec) {
			slice.push_back(v[i]);
		}
		return slice;
	}

	void storeOneWay_ET_EVt_VeT_VeVt_combis() {
		int count = 0;
		for (int T_id1 : T) {
			if (count % 1000 == 0) std::cout << count << " "; count++;
			Edge* e1 = rst->model->edgeVector[T_id1];
			if (!e1->silhouette) continue;
			Vertex* e1v1 = e1->vertices[0];
			Vertex* e1v2 = e1->vertices[1];
			std::vector<Vertex*> e1v = { e1v1, e1v2 };

			for (int T_id2 : T) {
				if (T_id2 == T_id1) continue;
				Edge* e2 = rst->model->edgeVector[T_id2];
				Vertex* e2v1 = e2->vertices[0];
				Vertex* e2v2 = e2->vertices[1];
				std::vector<Vertex*> e2v = { e2v1, e2v2 };

				bool e1_for_e2 = false;
				// e1_for e2v1, e2v2
				std::vector<bool> e1_for_e2v(2, false);
				// e1v1, e1v2 for e2
				std::vector<bool> e1v_for_e2(2, false);

				bool side[2] = { false, false };
				for (int i = 0; i < 2; i++) {
					e1_for_e2v[i] = e1->isSilhouetteForPos(e2v[i]->pos, side[i]);
					if (e1_for_e2v[i])e1_for_e2 = true;
				}
				if (!e1_for_e2) e1_for_e2 = (side[0] != side[1]);
				for (int i = 0; i < 2; i++) e1v_for_e2[i] = e1_for_e2 && e1v[i]->silhouette;
				if (!e1_for_e2) continue;

				std::vector<Ray> flat;
				for (Vertex* e1ver : e1v) {
					for (Vertex* e2ver : e2v) {
						if (e2ver->id != e1ver->id) flat.push_back(Ray(e1ver->pos, e2ver->pos));
					}
				}
				std::vector<glm::dvec3> pts;
				for (Ray& r : flat) pts.push_back(rst->window->rayIntersection(r));
				AxisAlignedPolygon poly(pts, rst->window->normal);
				if (!rst->window->boundingBoxIntersect(poly)) continue;
				ET.insert({ (uint32_t)e1->id, (uint32_t)e2->id });
				
				// remove silhouettes when vertices are equal
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						if (e1v[i] == e2v[j]) {
							e1_for_e2v[j] = false; e1v_for_e2[i] = false;
							e1_for_e2v[j] = false; e1v_for_e2[i] = false;
						}
					}
				}

				// rays
				// e1v1 e2v1, e1v1 e1v2
				// e1v2 e2v1, e1v2 e1v2
				std::vector<std::vector<Ray>> v_to_v(2, std::vector<Ray>(2));
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						v_to_v[i][j] = Ray(e1v[i]->pos, e2v[j]->pos);
					}
				}

				for (int i = 0; i < 2; i++) {
					if (e1_for_e2v[i] && !rst->window->intersectsPlaneFromLines(getCol(v_to_v, i)))
						e1_for_e2v[i] = false;
				}
				for (int i = 0; i < 2; i++) {
					if (e1v_for_e2[i] && !rst->window->intersectsPlaneFromLines(v_to_v[i]))
						e1v_for_e2[i] = false;
				}

				for (int i = 0; i < 2; i++) {
					if (e1_for_e2v[i]) EVt.insert({ (uint32_t)e1->id, (uint32_t)e2v[i]->id });
					if (e1v_for_e2[i]) VeT.insert({(uint32_t)e1v[i]->id, (uint32_t)e2->id});
					for (int j = 0; j < 2; j++) {
						if (e1_for_e2v[j] && e1v_for_e2[i]) VeVt.insert({ (uint32_t)e1v[i]->id, (uint32_t)e2v[j]->id });
					}
				}
			}
		}
	}

	void getTwoWayfromOneWay(std::unordered_set<std::vector<uint32_t>, hash>& oneWay1, std::unordered_set<std::vector<uint32_t>, hash>& oneWay2, std::unordered_set<std::vector<uint32_t>, hash>& twoWay, bool storeAsIs = false) {
		for (std::vector<uint32_t> ow1 : oneWay1) {
			if (oneWay2.count({ ow1[1], ow1[0] })) {
				if (ow1[0] < ow1[1] || storeAsIs) twoWay.insert(ow1);
				else twoWay.insert({ ow1[1], ow1[0] });
			}
		}
	}

	void getTwoWayfromOneWay(std::unordered_set<std::vector<uint32_t>, hash>& oneWay, std::unordered_set<std::vector<uint32_t>, hash>& twoWay) {
		for (std::vector<uint32_t> ow1 : oneWay) {
			if (oneWay.count({ ow1[1], ow1[0] })) {
				if (ow1[0] < ow1[1]) twoWay.insert(ow1);
				else twoWay.insert({ ow1[1], ow1[0] });
			}
		}
	}

	bool checkEEOcclusion(std::vector<uint32_t> ee) {
		Edge* e1 = rst->model->edgeVector[ee[0]];
		Edge* e2 = rst->model->edgeVector[ee[0]];
		std::vector<glm::dvec3> e1vec = e1->getVecPos();
		std::vector<glm::dvec3> e2vec = e2->getVecPos();
		Tetrahedron tetra(e1vec, e2vec, false);
		std::vector<std::set<int>> intersects;
		//int inter = 0;
		int r = 0;
		int raytracePrim;
		double raytraceDepth = 0;
		double check_offset = 1E-7;
		for (TriWedge& wedge : tetra.wedges) {
			for (Ray& ray : wedge.rays) {
				int inter = 0;
				double cast_offset = 0;
				double depth = std::max(ray.depthToIntersectionWithRay(e1->ray), ray.depthToIntersectionWithRay(e2->ray));
				if (glm::dot(ray.direction, rst->window->normal) > 0) ray.inverseDir();
				while (depth - cast_offset - check_offset < raytraceDepth) {
					ray.offsetByDepth(raytraceDepth + 0.001);
					cast_offset = raytraceDepth + 0.001;
					bool hit = rst->model->getIntersectionNoAcceleration(ray, raytracePrim, raytraceDepth);
					if (depth - cast_offset - check_offset > raytraceDepth) intersects[r].insert(raytracePrim);
				}
				if (intersects[r].empty()) return false;
				r++;
			}
		}
		for (int i1 : intersects[0]) {
			int countsame = 1;
			for (int i = 1; i < 4; i++) if (intersects[i].count(i1)) countsame++;
			if (countsame == 4) return true; // all lines intersect the same triangle

			Primitive* p1 = rst->model->triangles[i1];
			std::vector<Primitive*> neighbors = { p1 };
			for (int i = 1; i < 4; i++) {
				for (int i_inter : intersects[i]) {
					Primitive* p2 = rst->model->triangles[i_inter];
					for (Primitive* p : neighbors) {

					}
					if (p2->vertexNeighbor(p1)) { // and neighbors with others? or not?
						neighbors.push_back(p2);
						break;
					}
				}
			}
			if (neighbors.size() != 4) return false;
			// check if all have the same vertex
			std::set<Vertex*> sameverts = p1->sameVertices(neighbors[1]);
			       
			
			for (int i = 0; i < 4; i++) {
				int count = 0;
				// if all are edge neighbors with at least 2 others return true
				// else if not all are edge neigbours but other edges for shared vertex
				// are all convex, return true
			}

		}

		return false;
	}

	void storeEE_ET_EVe_EVt_VeT_VeVt_VeVe_combis() {
		storeOneWay_ET_EVt_VeT_VeVt_combis();

		std::cout << std::endl << "ET EE " << ET.size() << std::endl;
		getTwoWayfromOneWay(ET, EE);
		//std::vector<std::vector<uint32_t>> remove;
		//for (std::vector<uint32_t> ee : ET) {
		//	if (!rst->window->isInFrontOf(rst->model->edgeVector[ee[0]]->getVecPos(), rst->model->edgeVector[ee[1]]->getVecPos())) remove.push_back(ee);
		//}
		//for (std::vector<uint32_t> rem : remove) ET.erase(rem);

		std::cout << std::endl << "VeVt VeVe " << VeVt.size() << std::endl;
		getTwoWayfromOneWay(VeVt, VeVe);

		//remove = std::vector<std::vector<uint32_t>>();
		//for (std::vector<uint32_t> vv : VeVt) {
		//	if (!rst->window->isInFrontOf(rst->model->vertices[vv[0]]->pos, rst->model->vertices[vv[1]]->pos)) remove.push_back(vv);
		//}
		//for (std::vector<uint32_t> rem : remove) VeVt.erase(rem);

		std::cout << std::endl << "EVt VeT EVe " << EVt.size() + VeT.size() << std::endl;
		getTwoWayfromOneWay(EVt, VeT, EVe, true);
		//remove = std::vector<std::vector<uint32_t>>();
		//for (std::vector<uint32_t> ve : VeT) {
		//	if (!rst->window->isInFrontOf(rst->model->vertices[ve[0]]->pos, rst->model->edgeVector[ve[1]]->getVecPos())) remove.push_back(ve);
		//}
		//for (std::vector<uint32_t> rem : remove) VeT.erase(rem);

		//remove = std::vector<std::vector<uint32_t>>();
		//for (std::vector<uint32_t> ev : EVt) {
		//	if (!rst->window->isInFrontOf(rst->model->edgeVector[ev[0]]->getVecPos(), rst->model->vertices[ev[1]]->pos)) remove.push_back(ev);
		//}
		//for (std::vector<uint32_t> rem : remove) EVt.erase(rem);
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