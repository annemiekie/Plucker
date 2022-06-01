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
#include "extremalStabbingLine.h"
#include "splitSide.h"
#include "combiConfigurations.h"
#include <unordered_set>

class ESLFindAll {
public:
	std::vector<ESLCandidate> eslsNoVis;
	std::vector<ESLCandidate> esls;
	std::vector<Edge*> silhouetteEdges;
	bool checkHardCombis = false;
	std::vector<std::vector<ExtremalStabbingLine>> eslsPerPrimitive;
	std::vector<std::vector<Ray>> raysPerPrimitive;


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
	std::unordered_set<std::vector<uint32_t>, hash> EET;

	std::unordered_set<std::vector<uint32_t>, hash> VeVe;
	std::unordered_set<std::vector<uint32_t>, hash> VeVt;

	std::unordered_set<std::vector<uint32_t>, hash> EEVe;
	std::unordered_set<std::vector<uint32_t>, hash> EVeT;
	std::unordered_set<std::vector<uint32_t>, hash> EEVt;
	std::unordered_set<std::vector<uint32_t>, hash> EEET;
	std::unordered_set<std::vector<uint32_t>, hash> EEEE;

	std::unordered_set<std::vector<uint32_t>, hash> WT;
	std::unordered_set<std::vector<uint32_t>, hash> WE;
	std::unordered_set<std::vector<uint32_t>, hash> WEE;
	std::unordered_set<std::vector<uint32_t>, hash> WET;
	std::unordered_set<std::vector<uint32_t>, hash> WEET;
	std::unordered_set<std::vector<uint32_t>, hash> WEEE;
	std::unordered_set<std::vector<uint32_t>, hash> WVe;
	std::unordered_set<std::vector<uint32_t>, hash> WVt;
	std::unordered_set<std::vector<uint32_t>, hash> WEVe;
	std::unordered_set<std::vector<uint32_t>, hash> WEVt;
	std::unordered_set<std::vector<uint32_t>, hash> WVeT;

	std::unordered_set<std::vector<uint32_t>, hash> VwVe;
	std::unordered_set<std::vector<uint32_t>, hash> VwVt;
	std::unordered_set<std::vector<uint32_t>, hash> VwE;
	std::unordered_set<std::vector<uint32_t>, hash> VwT;
	std::unordered_set<std::vector<uint32_t>, hash> VwEE;
	std::unordered_set<std::vector<uint32_t>, hash> VwET;


	// all vertices from viable edges
	//std::vector<std::vector<uint32_t>> VtTT; // in plane through vertex

	//std::vector<std::vector<uint32_t>> P; // plane of primitive (one per facade)
	//std::vector<std::vector<uint32_t>> EP;
	//std::vector<std::vector<uint32_t>> EEP;

	//std::vector<std::vector<uint32_t>> TT; // 

	CombiConfigurations combiS;

	bool findAllEsl() {
		markPotentialSilhouettes();
		markSplitLines();

		storeT();											// SSST SSSE
		storeVt();											// SSVt SSVe
		std::cout << "EE_ET_EVe_EVt_VeT_VeVt_VeVe " << T.size() << "x" << T.size() << std::endl;
		storeEE_ET_EVe_EVt_VeT_VeVt_VeVe_combis();			// SSEE SSET SEVt SVeT SEVe VeVe VeVt
		//std::cout << "EEE_EET " << EE.size() << "x" << T.size() << std::endl;
		//storeEEE_EET_combis();								// SEEE SEET
		//std::cout << "EEVe_EVeT " << EVe.size() << "x" << T.size() << std::endl;
		//storeEEVe_EVeT_combis();							// EEVe EVeT
		//std::cout << "EEEE_EEET " << EEE.size() << "x" << T.size() << std::endl;
		//storeEEEE_EEET_combis();							// EEEE EEET
		std::cout << "WE_WT_WVt_WVe " << T.size() << "x" << 4 << std::endl;
		storeWE_WT_WVt_WVe_combis();						// WSSE WSST WSVe WSVt
		std::cout << "VwVe_VwVt_VwE_VwT " << T.size() << "x" << 4 << std::endl;
		storeVwVe_VwVt_VwE_VwTcombis();						// VwVe VwVt VwST VwSE
		std::cout << "VwEE_VwET " << VwE.size() << "x" << T.size() << std::endl;
		storeVwEE_VwETcombis();								// VwEE VwET
		std::cout << "WEVt_WEVe " << WE.size() << "x" << Vt.size() << std::endl;
		storeWEVt_WEVe_combis();							// WEVt WEve
		//std::cout << "WVeT " << WVe.size() << "x" << T.size() << std::endl;
		//storeWVeTcombis();									// WVeT
		//std::cout << "WEE_WET " << WE.size() << "x" << T.size() << std::endl;
		//storeWEE_WET_combis();								// WSEE	WSET
		//std::cout << "WEEE_WEET " << WEE.size() << "x" << T.size() << std::endl;
		//storeWEEE_WEET_combis();							// WEEE WEET		

		// check if inplane stuff is covered! or add seperately
		// or make sure that if edgeA is in plane of 1 of 2 silh triangles of edgeB, edgeB is silh for edgeA
		// especially extra check to see if part of same primitive -> definitely in plane
		makeESLsFromCombi();
		return true;
	}

	Line4 createEEVLine(Edge* e1, Edge* e2, glm::dvec3 v) {
		TriWedge e1v(v, e1->vertices[0]->pos, e1->vertices[1]->pos, glm::vec3(0, 0, 0));
		glm::dvec3 to = e1v.rayIntersection(e2->ray);
		return Line4(v, to);
	}

	Line4 createEEVLine(std::vector<glm::dvec3> e1, Edge* e2, Vertex* v) {
		TriWedge e1v(v->pos, e1[0], e1[1], glm::vec3(0, 0, 0));
		glm::dvec3 from = e1v.rayIntersection(e2->ray);
		Line4 ray(from, v->pos);
		if (glm::dot(ray.direction, rst->window->normal) > 0) ray.inverseDir();
		double depth = rst->window->rayIntersectionDepth(ray);
		ray.offsetByDepth(depth);
		return ray;
	}

	Line4 createEEVLine(Edge* e1, Edge* e2, Vertex* v) {
		TriWedge e1v(v->pos, e1->vertices[0]->pos, e1->vertices[1]->pos, glm::vec3(0, 0, 0));
		glm::dvec3 from = v->pos;
		glm::dvec3 to = e1v.rayIntersection(e2->ray);
		Line4 ray(from, to);
		if (glm::dot(ray.direction, rst->window->normal) > 0) ray.inverseDir();
		double depth = rst->window->rayIntersectionDepth(ray);
		ray.offsetByDepth(depth);
		return ray;
	}

	Ray createVVLine(Vertex* v1, Vertex* v2) {
		return Ray(v1->pos, v2->pos);
	}

	Ray createVVLine(glm::dvec3 v1, Vertex* v2) {
		return Ray(v1, v2->pos);
	}

	std::vector<Line4> createEEEELine(std::vector<Edge*> lines) {
		std::vector<Ray> rays;
		for (Edge* e : lines) rays.push_back(e->ray);
		return Lines4Finder::find(rays);
	}

	bool vertexSilhouetteForRay(Vertex* v, Ray& ray) {
		int count = 0;
		for (Edge* e : v->edges) {
			if (e->isSilhouetteForRay(ray)) count++;
			if (count == 2) return true;
		}
		return false;
	}

	void makeESLsFromCombi() {
		raysPerPrimitive.resize(rst->model->primsize);
		// first everything without split lines
		// VwVt
		//std::cout << "make VwVt " << VwVt.size() << std::endl;
		//for (std::vector<uint32_t> vwvt : VwVt) {
		//	glm::dvec3 vw = rst->window->vertices[vwvt[0]];
		//	Vertex* vt = rst->model->vertices[vwvt[1]];
		//	Line4 esl = createVVLine(vw, vt);
		//	// window check done
		//	// raytrace check
		//	double depth = esl.depthToPointOnRay(vt->pos);
		//	double embreedepth; int embreeprim;
		//	bool hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	//if (!hit) continue;
		//	if (depth - 1E-7 > embreedepth) continue;
		//	for (Primitive* p : vt->triangles) {
		//		if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
		//	}
		//}

		//// VwVe
		//std::cout << "make VwVe " << VwVe.size() << std::endl;
		//for (std::vector<uint32_t> vwve : VwVe) {
		//	glm::dvec3 vw = rst->window->vertices[vwve[0]];
		//	Vertex* ve = rst->model->vertices[vwve[1]];
		//	Ray esl = createVVLine(vw, ve);
		//	// window check done
		//	// raytrace check
		//	double depth = esl.depthToPointOnRay(ve->pos);
		//	double embreedepth; int embreeprim;
		//	bool hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	//rst->model->getIntersectionEmbree(esl, embreeprim, embreedepth, true);
		//	if (depth - 1E-7 > embreedepth)  continue;
		//	for (Primitive* p : ve->triangles) {
		//		if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
		//	}
		//	if (!vertexSilhouetteForRay(ve, esl) || !hit) continue;
		//	esl.offsetByDepth(depth + 0.001);
		//	rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);// , true);
		//	if (embreeprim >= 0) raysPerPrimitive[embreeprim].push_back(esl);
		//}


		//// VwEE
		//std::cout << "make VwEE " << VwEE.size() << std::endl;
		//for (std::vector<uint32_t> vwee : VwEE) {
		//	glm::dvec3 vw = rst->window->vertices[vwee[0]];
		//	Edge* e1 = rst->model->edgeVector[vwee[1]];
		//	Edge* e2 = rst->model->edgeVector[vwee[2]];
		//	std::vector<Edge*> edges = { e1, e2 };

		//	Line4 esl = createEEVLine(e1, e2, vw);
		//	// window check done
		//	// intersection check
		//	//if (!e1->intersectsRay(esl) || !e2->intersectsRay(esl)) continue;
		//	bool e1infrontofe2 = rst->window->isInFrontOf(esl.pointOfintersectWithRay(e1->ray), esl.pointOfintersectWithRay(e2->ray));
		//	uint16_t elast = e1infrontofe2 ? 1 : 0;

		//	// silhouette check
		//	bool e1silh = e1->isSilhouetteForRay(esl);
		//	bool e2silh = e2->isSilhouetteForRay(esl);
		//	if (!e1silh && !e2silh) continue;
		//	if (!e1silh && elast == 1) continue;
		//	if (!e2silh && elast == 0) continue;
		//	if (!e1silh || !e2silh) {
		//		VwET.insert({ vwee[0], (uint32_t)edges[1-elast]->id , (uint32_t)edges[elast]->id });
		//		continue;
		//	}

		//	// raytrace check
		//	double depthEfirst = esl.depthToIntersectionWithRay(edges[1-elast]->ray);
		//	double depthElast = esl.depthToIntersectionWithRay(edges[elast]->ray);
		//	double embreedepth; int embreeprim;
		//	bool hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	if (depthEfirst - 1E-7 > embreedepth) continue;
		//	esl.offsetByDepth(depthEfirst + 0.001);
		//	depthElast -= (depthEfirst + 0.001);
		//	hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	if (depthElast - 1E-7 > embreedepth) continue;

		//	// add to primitives of last edge
		//	for (Primitive* p : edges[elast]->triangles) {
		//		if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
		//	}
		//	if (!hit) continue;
		//	esl.offsetByDepth(depthElast + 0.001);
		//	rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);// , true);
		//	if (embreeprim >= 0) raysPerPrimitive[embreeprim].push_back(esl);
		//}

		//// VwET
		//std::cout << "make VwET " << VwET.size() << std::endl;
		//for (std::vector<uint32_t> vwet : VwET) {
		//	glm::dvec3 vw = rst->window->vertices[vwet[0]];
		//	Edge* e = rst->model->edgeVector[vwet[1]];
		//	Edge* t = rst->model->edgeVector[vwet[2]];

		//	Line4 esl = createEEVLine(e, t, vw);
		//	// window check done
		//	// silhouette check
		//	if (!e->isSilhouetteForRay(esl)) continue;
		//	// intersection check (not needed?)
		//	// raytrace check
		//	double depthE = esl.depthToIntersectionWithRay(e->ray);
		//	double depthT = esl.depthToIntersectionWithRay(t->ray);
		//	if (depthT < depthE) continue;
		//	double embreedepth; int embreeprim;
		//	bool hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	if (depthE - 1E-7 > embreedepth) continue;
		//	esl.offsetByDepth(depthE + 0.001);
		//	depthT -= (depthE + 0.001);
		//	hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	if (depthT - 1E-7 > embreedepth) continue;

		//	for (Primitive* p : t->triangles) {
		//		if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
		//	}
		//}

		//// WEVt
		//std::cout << "make WEVt " << WEVt.size() << std::endl;
		//for (std::vector<uint32_t> wevt : WEVt) {
		//	std::vector<glm::dvec3> w = { rst->window->vertices[wevt[0]],  rst->window->vertices[(wevt[0] + 1)%rst->window->vertices.size()] };
		//	Edge* e = rst->model->edgeVector[wevt[1]];
		//	Vertex* vt = rst->model->vertices[wevt[2]];

		//	Line4 esl = createEEVLine(w, e, vt);
		//	// window check done
		//	// silhouette check
		//	if (!e->isSilhouetteForRay(esl)) continue;
		//	// intersection check (not needed?)
		//	// raytrace check
		//	double depthE = esl.depthToIntersectionWithRay(e->ray);
		//	double depthVt = esl.depthToPointOnRay(vt->pos);
		//	if (depthE > depthVt) continue;
		//	double embreedepth; int embreeprim;
		//	bool hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	if (depthE - 1E-7 > embreedepth) continue;

		//	esl.offsetByDepth(depthE + 0.001);
		//	depthVt -= (depthE + 0.001);
		//	hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
		//	if (depthVt - 1E-7 > embreedepth) continue;

		//	for (Primitive* p : vt->triangles) {
		//		if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
		//	}
		//}

		// WEVe
		std::cout << "make WEVe " << WEVe.size() << std::endl;
		for (std::vector<uint32_t> weve : WEVe) {
			std::vector<glm::dvec3> w = { rst->window->vertices[weve[0]],  rst->window->vertices[(weve[0] + 1) % rst->window->vertices.size()] };
			Edge* e = rst->model->edgeVector[weve[1]];
			Vertex* ve = rst->model->vertices[weve[2]];

			Line4 esl = createEEVLine(w, e, ve);
			// window check done
			// silhouette check
			if (!e->isSilhouetteForRay(esl) || !vertexSilhouetteForRay(ve, esl)) continue;
			// intersection check (not needed?)
			// raytrace check
			double depthE = esl.depthToIntersectionWithRay(e->ray);
			double depthVe = esl.depthToPointOnRay(ve->pos);
			std::vector<double> depth;
			if (depthE < depthVe) depth = { depthE, depthVe };
			else depth = { depthVe, depthE };

			double embreedepth; int embreeprim;
			bool hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
			if (depth[0] - 1E-7 > embreedepth) continue;

			esl.offsetByDepth(depth[0] + 0.001);
			depth[1] -= (depth[0] + 0.001);
			hit = rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);
			if (depth[1] - 1E-7 > embreedepth) continue;

			if (depthE < depthVe) {
				for (Primitive* p : ve->triangles) {
					if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
				}
			}
			else {
				for (Primitive* p : e->triangles) {
					if (glm::dot(p->normal, esl.direction) < 0) raysPerPrimitive[p->id].push_back(esl);
				}
			}
			if (!hit) continue;
			esl.offsetByDepth(depth[1] + 0.001);
			rst->model->getIntersectionNoAcceleration(esl, embreeprim, embreedepth);// , true);
			if (embreeprim >= 0) raysPerPrimitive[embreeprim].push_back(esl);
		}
		// WEET

		// WEEE

		// EEVe

		// EEEE

		// EEET

		// EEVt

		// VeVt

		// VeVe

	}

	bool markSplitLines() {
		for (Split split : rst->splitters) {
			split.edge->splitline = true;
			for (Vertex* v : split.edge->vertices) v->splitline = true;
		}
		return true;
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
			if (count > 0 && count < rst->window->vertices.size()) e->silhouette = true;
			if (e->silhouette) {
				for (Vertex* v : e->vertices) v->silhouette = true;
			}
		}
		return true;
	}

	void storeVt() {
		for (Vertex* v : rst->model->vertices) {
			for (Primitive* prim : v->triangles) {
				if (glm::dot(rst->window->normal, prim->normal) > 0 ||
					rst->window->intersectsPlaneFromLines(prim->getRayVector())) {
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
					rst->window->intersectsPlaneFromLines(prim->getRayVector())) {
					T.insert(e->id);
					break;
				}
			}
		}
	}

	void storeWVeTcombis() {
		for (std::vector<uint32_t> vec : WVe) {
			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			std::vector<glm::dvec3> w1v = { rst->window->vertices[id1], rst->window->vertices[(id1 + 1) % rst->window->vertices.size()] };
			Vertex* v2 = rst->model->vertices[id2];
			TriWedge w1v2(v2->pos, w1v[0], w1v[1], glm::vec3(0, 0, 0));

			for (uint32_t id3 : T) {
				Edge* e3 = rst->model->edgeVector[id3];
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				glm::dvec3 pt;
				if (!w1v2.raySegmentIntersection(e3v[0], e3v[1], pt)) continue;
				if (findVecInSet(WEVe, { id1, id3, id2 })) continue;
				bool W1T3 = findVecInSet(WT, {id1, id3});
				bool Ve2T3 = findVecInSet(EVe, {id2, id3});
				if (W1T3 && Ve2T3)  WVeT.insert({ id1, id2, id3 });
				bool W1E3 = findVecInSet(WE, { id1, id3 });
				bool E3Ve2 = findVecInSet(EVe, { id2, id3 });

				if (rst->window->isInFrontOf(v2->pos, e3v) && (W1T3 || W1E3) && (Ve2T3 || E3Ve2)) WVeT.insert({ id1, id2, id3 });
			}
		}
	}

	void storeWEVt_WEVe_combis() {
		for (std::vector<uint32_t> vec : WE) {
			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e2 = rst->model->edgeVector[id2];

			std::vector<glm::dvec3> w1v = { rst->window->vertices[id1], rst->window->vertices[(id1 + 1) % rst->window->vertices.size()] };
			std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };

			for (uint32_t id3 : Vt) {
				Vertex* v3 = rst->model->vertices[id3];

				TriWedge e1v3(v3->pos, w1v[0], w1v[1], glm::vec3(0, 0, 0));
				glm::dvec3 pt;
				if (!e1v3.raySegmentIntersection(e2v[0], e2v[1], pt)) continue;

				bool W1Ve3 = false;
				bool E2Ve3 = false;
				if (v3->silhouette) {
					W1Ve3 = findVecInSet(WVe, std::vector<uint32_t>{id1, id3});
					E2Ve3 = findVecInSet(EVe, std::vector<uint32_t>{id2, id3});
					if (W1Ve3 && E2Ve3) {
						WEVe.insert({ id1, id2, id3 });
						continue;
					}
				}
				if (!rst->window->isInFrontOf(e2v, v3->pos)) continue;
				bool W1Vt3 = findVecInSet(WVt, std::vector<uint32_t>{id1, id3});
				bool E2Vt3 = findVecInSet(EVt, std::vector<uint32_t>{id2, id3});
				if ((W1Vt3 || W1Ve3) && (E2Ve3 || E2Vt3)) WEVt.insert({ id1, id2, id3 });
			}
		}
	}

	void storeVwEE_VwETcombis() {
		for (std::vector<uint32_t> vec : VwE) {
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
					E2E3 = findVecInSet(EE, { id2, id3 });
					Vw1E3 = findVecInSet(VwE, { id1, id3 });
					if (E2E3 && Vw1E3) {
						VwEE.insert({id1, id2, id3});
						continue;
					}
				}
				if (!rst->window->isInFrontOf(e2v, e3v)) continue;
				bool E2T3 = false;
				bool Vw1T3 = false;
				E2T3 = findVecInSet(ET, { id2, id3 });
				Vw1T3 = findVecInSet(VwT, { id1, id3 });
				if ((E2T3 || E2E3) && (Vw1E3 || Vw1T3)) VwET.insert({ id1, id2, id3 });
			}
		}
	}

	bool findVecInSet(std::unordered_set<std::vector<uint32_t>, hash>& set, std::vector<uint32_t> vec) {
		return set.count(vec);
		//for (std::vector<uint32_t> v : set) {
		//	if (v == vec) return true;
		//}
		//return false;
	}

	void storeWEE_WET_combis() {
		for (std::vector<uint32_t> vec : WE) {
			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e2 = rst->model->edgeVector[id2];

			std::vector<glm::dvec3> w1v = { rst->window->vertices[id1], rst->window->vertices[(id1+1) % rst->window->vertices.size()] };
			std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };
			Tetrahedron unboundTetra(w1v, e2v);

			for (uint32_t id3 : T) {
				Edge* e3 = rst->model->edgeVector[id3];
				std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
				if (!unboundTetra.segmInTetra(e3v[0], e3v[1], 1E-7)) continue;

				bool W1T3 = false;
				bool E2T3 = false;
				bool W1E3 = false;
				bool E2E3 = false;

				if (e3->silhouette && id3 > id2) {
					E2E3 = findVecInSet(EE, std::vector<uint32_t>{id2, id3});
					W1E3 = findVecInSet(WE, std::vector<uint32_t>{id1, id3});
					if (W1E3 && E2E3) {
						WEE.insert({ id1, id2, id3 });
						continue;
					}
				}
				if (!rst->window->isInFrontOf(e2v, e3v)) continue;
				W1T3 = findVecInSet(WT, std::vector<uint32_t>{id1, id3});
				E2T3 = findVecInSet(ET, std::vector<uint32_t>{id2, id3});
				if ((W1T3 && E2T3) || (W1T3 && E2E3) || (W1E3 && E2T3)) WET.insert({ id1, id2, id3 });
			}
		}
	}

	void storeWEEE_WEET_combis() {
		for (std::vector<uint32_t> vec : WEE) {
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
					W1E2E4 = findVecInSet(WEE, std::vector<uint32_t>{id1, id2, id4});
					W1E3E4 = findVecInSet(WEE, std::vector<uint32_t>{id1, id3, id4});
					E2E3E4 = findVecInSet(EEE, std::vector<uint32_t>{id2, id3, id4});
					if (W1E2E4 && W1E3E4 && E2E3E4) {
						WEEE.insert({ id1, id2, id3, id4 });
						continue;
					}
				}
				if (!rst->window->isInFrontOf(e2->getVecPos(), e4->getVecPos()) || !rst->window->isInFrontOf(e3->getVecPos(), e4->getVecPos())) continue;
				bool W1E2T4 = findVecInSet(WT, std::vector<uint32_t>{id1, id2, id4});
				bool W1E3T4 = findVecInSet(ET, std::vector<uint32_t>{id1, id3, id4});
				bool E2E3T4 = findVecInSet(EET, std::vector<uint32_t>{id2, id3, id4});
				if ((E2E3T4 || E2E3E4) && (W1E2E4 || W1E2T4) && (W1E3E4 || W1E3T4)) WEET.insert({ id1, id2, id3 });
			}
		}
	}

	void storeVwVe_VwVt_VwE_VwTcombis() {
		for (uint32_t id1 : Vt) {
			Vertex* v2 = rst->model->vertices[id1];
			for (uint32_t i = 0; i < rst->window->vertices.size(); i++) {
				bool W1aVe2 = false;
				bool W1bVe2 = false;
				uint32_t j = (i + 1) % rst->window->vertices.size();
				if (v2->silhouette) {
					W1aVe2 = findVecInSet(WVe, { i, id1 });
					W1bVe2 = findVecInSet(WVe, { j, id1 });
					if (W1aVe2 && W1bVe2) {
						VwVe.insert({ i, id1 });
						continue;
					}
				}
				bool W1aVt2 = false;
				bool W1bVt2 = false;
				W1aVt2 = findVecInSet(WVt, { i, id1 });
				W1bVt2 = findVecInSet(WVt, { j, id1 });
				if ((W1aVe2 || W1aVt2) && (W1bVe2 || W1bVt2)) VwVt.insert({ i, id1 });

			}
		}
		for (uint32_t id1 : T) {
			Edge* e = rst->model->edgeVector[id1];
			for (uint32_t i = 0; i < rst->window->vertices.size(); i++) {
				bool W1aE2 = false;
				bool W1bE2 = false;
				uint32_t j = (i + 1) % rst->window->vertices.size();
				if (e->silhouette) {
					W1aE2 = findVecInSet(WE, { i, id1 });
					W1bE2 = findVecInSet(WE, { j, id1 });
					if (W1aE2 && W1bE2) {
						VwE.insert({ i, id1 });
						continue;
					}
				}
				bool W1aT2 = false;
				bool W1bT2 = false;
				W1aT2 = findVecInSet(WT, { i, id1 });
				W1bT2 = findVecInSet(WT, { j, id1 });
				if ((W1aE2 || W1aT2) && (W1bE2 || W1bT2)) VwT.insert({ i, id1 });
			}
		}
	}

	void storeWE_WT_WVt_WVe_combis() {
		for (uint32_t id1 : T) {
			Edge* e1 = rst->model->edgeVector[id1];

			for (uint32_t i = 0; i < rst->window->vertices.size(); i++) {
				glm::dvec3 v1 = rst->window->vertices[i];
				glm::dvec3 v2 = rst->window->vertices[(i+1)%rst->window->vertices.size()];

				if (e1->silhouette && e1->isSilhouetteForPrim({ v1, v2 })) {
					WE.insert({ i, id1 });
					WVe.insert({i, uint32_t(e1->vertices[0]->id)});
					WVe.insert({ i, uint32_t(e1->vertices[1]->id) });
				}
				else {
					WT.insert({ i, id1 });
					WVt.insert({ i, uint32_t(e1->vertices[0]->id) });
					WVt.insert({ i, uint32_t(e1->vertices[1]->id) });
				}
			}
		}
		std::vector<std::vector<uint32_t>> remove;
		for (std::vector<uint32_t> WVtvec : WVt) {
			if (findVecInSet(WVe, WVtvec)) remove.push_back(WVtvec);
		}
		for (std::vector<uint32_t> rem : remove) WVt.erase(rem);
	}

	void storeEEEE_EEET_combis() {
		for (std::vector<uint32_t> vec : EEE) {
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
					E1E2E4 = findVecInSet(EEE, { id1, id2, id4 });
					E2E3E4 = findVecInSet(EEE, { id2, id3, id4 });
					E1E3E4 = findVecInSet(EEE, { id1, id3, id4 });
					if (E1E2E4 && E2E3E4 && E1E3E4) {
						EEEE.insert({ id1, id2, id3, id4 });
						continue;
					}
				}

				if (!(rst->window->isInFrontOf(e1->getVecPos(), e4->getVecPos()) && rst->window->isInFrontOf(e2->getVecPos(), e4->getVecPos()) && rst->window->isInFrontOf(e3->getVecPos(), e4->getVecPos()))) continue;

				bool E1E2T4 = findVecInSet(EET, { id1, id2, id4 });
				if (!E1E2T4 && !E1E2E4) continue;
				bool E2E3T4 = findVecInSet(EET, { id2, id3, id4 });
				if (!E2E3T4 && !E2E3E4) continue;
				bool E1E3T4 = findVecInSet(EET, { id1, id3, id4 });
				if (!E1E3E4 && !E1E3T4) continue;
				EEET.insert({ id1, id2, id3, id4 });
			}
		}
	}

	void storeEEVt_combis() {
		for (std::vector<uint32_t> vec : EE) {
			uint32_t id1 = vec[0];
			uint32_t id2 = vec[1];
			Edge* e1 = rst->model->edgeVector[vec[0]];
			Edge* e2 = rst->model->edgeVector[vec[1]];
			std::vector<glm::dvec3> e1v = e1->getVecPos();
			std::vector<glm::dvec3> e2v = e2->getVecPos();

			for (uint32_t id3 : Vt) {
				Vertex* v3 = rst->model->vertices[id3];

				if (e1->vertices[0] == v3 || e1->vertices[1] == v3 ||
					e2->vertices[0] == v3 || e2->vertices[1]) continue;

				TriWedge e1v3(v3->pos, e1v[0], e1v[1], glm::vec3(0, 0, 0));
				glm::dvec3 pt;
				if (!e1v3.raySegmentIntersection(e2v[0], e2v[1], pt)) continue;

				if (findVecInSet(EEVe, { id1, id2, id3 })) continue;

				bool E1Vt3 = findVecInSet(EVt, { id1, id3 });
				bool E2Vt3 = findVecInSet(EVt, { id2, id3 });
				if (E1Vt3 && E2Vt3) {
					EEVt.insert({ id1, id2, id3 });
					continue;
				}
				bool E1Ve3 = findVecInSet(EVe, { id1, id3 });
				bool E2Ve3 = findVecInSet(EVe, { id2, id3 });
				if (!rst->window->isInFrontOf(e1v, v3->pos) || !rst->window->isInFrontOf(e2v, v3->pos)) continue;
				if ((E1Vt3 || E1Ve3) && (E2Vt3 || E2Ve3)) EEVt.insert({ id1, id2, id3 });
			}
		}
	}

	void storeEEVe_EVeT_combis() {
		// check for shared vertices??
		for (std::vector<uint32_t> vec : EVe) {
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

				bool E1T3 = false;
				bool V2T3 = false;
				bool E1E3 = false;
				bool V2E3 = false;

				if (e3->silhouette && id3 > id1) {
					E1E3 = findVecInSet(EE, { id1, id3 });
					V2E3 = findVecInSet(EVe, { id3, id2 });
					if (E1E3 && V2E3) {
						EEVe.insert({ id1, id2, id3 });
						continue;
					}
				}

				if (rst->window->isInFrontOf(e1v, e3v) && rst->window->isInFrontOf(v2->pos, e3v)) {
					E1T3 = findVecInSet(ET, { id1, id3 });
					V2T3 = findVecInSet(VeT, { id2, id3 });
					if ((E1T3 || E1E3) && (V2E3 || V2T3)) EVeT.insert({ id1, id2, id3 });					
				}
			}
		}
	}

	void storeEEE_EET_combis() {

		for (std::vector<uint32_t> vec : EE) {
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

				std::vector<uint32_t> E1E3vec{ id1, id3 };
				std::vector<uint32_t> E2E3vec{ id2, id3 };

				bool E1T3 = false;
				bool E2T3 = false;
				bool E1E3 = false;
				bool E2E3 = false;

				if (e3->silhouette && id3 > id2) {
					for (std::vector<uint32_t> vec2 : EE) {
						if (vec2 == E1E3vec) E1E3 = true;
						if (vec2 == E2E3vec) E2E3 = true;
						if (E1E3 && E2E3) break;
					}
					if (E1E3 && E2E3) {
						EEE.insert({ id1, id2, id3 });
						continue;
					}
				}

				if (!(rst->window->isInFrontOf(e1v, e3v) && rst->window->isInFrontOf(e2v, e3v))) continue;
				for (std::vector<uint32_t> vec2 : ET) {
					if (vec2 == std::vector<uint32_t>{id1, id3}) E1T3 = true;
					else if (vec2 == std::vector<uint32_t>{id2, id3}) E2T3 = true;
					if ((E1E3 && E2T3) || (E2E3 && E1T3) || (E1T3 && E2T3)) break;
				}
				if ((E1T3 || E1E3) && (E2T3 || E2E3)) EET.insert({ id1, id2, id3 });
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
		for (int T_id1 : T) {
			Edge* e1 = rst->model->edgeVector[T_id1];
			Vertex* e1v1 = e1->vertices[0];
			Vertex* e1v2 = e1->vertices[1];
			std::vector<Vertex*> e1v = { e1v1, e1v2 };

			for (int T_id2 : T) {
				if (T_id2 == T_id1) continue;
				Edge* e2 = rst->model->edgeVector[T_id2];
				Vertex* e2v1 = e2->vertices[0];
				Vertex* e2v2 = e2->vertices[1];
				std::vector<Vertex*> e2v = { e2v1, e2v2 };
				std::vector<std::vector <Vertex* >> v = { e1v, e2v };

				bool e1_for_e2 = false;
				// e1_for e2v1, e2v2
				std::vector<bool> e1_for_e2v(2, false);

				// e1v1, e1v2 for e2
				std::vector<bool> e1v_for_e2(2, false);

				if (e1->silhouette) {
					bool side[2] = { false, false };
					for (int i = 0; i < 2; i++) {
						e1_for_e2v[i] = e1->isSilhouetteForPos(v[1][i]->pos, side[i]);
						if (e1_for_e2v[i])e1_for_e2 = true;
					}
					if (!e1_for_e2) e1_for_e2 = (side[0] != side[1]);
					for (int i = 0; i < 2; i++) e1v_for_e2[i] = e1_for_e2 && e1v[i]->silhouette;
				}
				if (!e1_for_e2) continue;
				ET.insert({ (uint32_t)e1->id, (uint32_t)e2->id });

				// remove silhouettes when vertices are equal
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						if (v[0][i] == v[1][j]) {
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
						v_to_v[i][j] = Ray(v[0][i]->pos, v[1][j]->pos);
					}
				}

				// remove everything not intersecting window
				// 
				//fix this -->if (!rst->window->intersectsPlaneFromLines(flatten(v_to_v))) continue;
				for (int i = 0; i < 2; i++) {
					if (e1_for_e2v[i] && !rst->window->intersectsPlaneFromLines(getCol(v_to_v, i)))
						e1_for_e2v[i] = false;
				}
				for (int i = 0; i < 2; i++) {
					if (e1v_for_e2[i] && !rst->window->intersectsPlaneFromLines(v_to_v[i]))
						e1v_for_e2[i] = false;
				}

				for (int i = 0; i < 2; i++) {
					if (e1_for_e2v[i]) EVt.insert({ (uint32_t)e1->id, (uint32_t)v[1][i]->id });
					if (e1v_for_e2[i]) VeT.insert({(uint32_t)v[0][i]->id, (uint32_t)e2->id});
					for (int j = 0; j < 2; j++) {
						if (e1_for_e2v[j] && e1v_for_e2[i]) VeVt.insert({ (uint32_t)v[0][i]->id, (uint32_t)v[1][j]->id });
					}
				}
			}
		}
	}

	void getTwoWayfromOneWay(std::unordered_set<std::vector<uint32_t>, hash>& oneWay1, std::unordered_set<std::vector<uint32_t>, hash>& oneWay2, std::unordered_set<std::vector<uint32_t>, hash>& twoWay) {
		for (std::vector<uint32_t> ow1 : oneWay1) {
			for (std::vector<uint32_t> ow2 : oneWay2) {
				if (ow1[0] == ow2[1] && ow1[1] == ow2[0]) twoWay.insert(ow1);
			}
		}
		for (std::vector<uint32_t> tw : twoWay) {
			oneWay1.erase(tw);
			oneWay2.erase({ tw[1], tw[0] });
		}
	}

	void getTwoWayfromOneWay(std::unordered_set<std::vector<uint32_t>, hash>& oneWay, std::unordered_set<std::vector<uint32_t>, hash>& twoWay) {
		for (std::vector<uint32_t> ow1 : oneWay) {
			for (std::vector<uint32_t> ow2 : oneWay) {
				if (ow2[0] < ow1[0]) continue;
				if (ow1[0] == ow2[1] && ow1[1] == ow2[0]) twoWay.insert(ow1);
			}
		}
		for (std::vector<uint32_t> tw : twoWay) {
			oneWay.erase({ tw[0], tw[1] });
			oneWay.erase({ tw[1], tw[0] });
		}
	}

	void storeEE_ET_EVe_EVt_VeT_VeVt_VeVe_combis() {
		storeOneWay_ET_EVt_VeT_VeVt_combis();

		getTwoWayfromOneWay(ET, EE);
		std::vector<std::vector<uint32_t>> remove;
		for (std::vector<uint32_t> ee : ET) {
			if (!rst->window->isInFrontOf(rst->model->edgeVector[ee[0]]->getVecPos(), rst->model->edgeVector[ee[1]]->getVecPos())) remove.push_back(ee);
		}
		for (std::vector<uint32_t> rem : remove) ET.erase(rem);

		getTwoWayfromOneWay(VeVt, VeVe);
		remove = std::vector<std::vector<uint32_t>>();
		for (std::vector<uint32_t> vv : VeVt) {
			if (!rst->window->isInFrontOf(rst->model->vertices[vv[0]]->pos, rst->model->vertices[vv[1]]->pos)) remove.push_back(vv);
		}
		for (std::vector<uint32_t> rem : remove) VeVt.erase(rem);

		getTwoWayfromOneWay(EVt, VeT, EVe);
		remove = std::vector<std::vector<uint32_t>>();
		for (std::vector<uint32_t> ve : VeT) {
			if (!rst->window->isInFrontOf(rst->model->vertices[ve[0]]->pos, rst->model->edgeVector[ve[1]]->getVecPos())) remove.push_back(ve);
		}
		for (std::vector<uint32_t> rem : remove) VeT.erase(rem);

		remove = std::vector<std::vector<uint32_t>>();
		for (std::vector<uint32_t> ev : EVt) {
			if (!rst->window->isInFrontOf(rst->model->edgeVector[ev[0]]->getVecPos(), rst->model->vertices[ev[1]]->pos)) remove.push_back(ev);
		}
		for (std::vector<uint32_t> rem : remove) EVt.erase(rem);
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