#pragma once
#pragma once
#include "lineThroughFourLines.h"
#include "edge.h"
#include "vertex.h"
#include "splitSide.h"
#include "triWedge.h"
#include "axisAlignedPolygon.h"
#include "rst.h"
#include "findLineThroughFour.h"

struct ESLCandidatePre {
	std::vector<Line4> rays;
	int nrVertices = 0;
	//bool allSilhouette = false;

	RaySpaceTree* rst;
	glm::dvec3 Vw;
	std::vector<glm::dvec3> W;
	std::vector<Vertex*> Ve;
	std::vector<Edge*> E;
	//Vertex* Vt;
	//Edge* T;
	std::vector<uint32_t> indices;
	int nrVw = 0; int nrW = 0; int nrVe = 0; int nrE = 0; //int nrVt = 0; int nrT = 0;

	std::vector<double> depth;
	std::vector<std::set<int>> intprims;

	int id = -1;

	ESLCandidatePre(std::vector<uint32_t>& indices, int nrVertices, //bool allsilh, 
					RaySpaceTree* rst,
					int nrVw, int nrW, int nrVe, int nrE) ://, int nrVt, int nrT) :
					indices(indices), nrVertices(nrVertices), //allSilhouette(allsilh), 
					rst(rst),
					nrVw(nrVw), nrW(nrW), nrE(nrE), nrVe(nrVe) { //, nrVt(nrVt), nrT(nrT) {
		int c = 0;
		if (nrVw) Vw = rst->window->vertices[indices[c++]];
		else if (nrW) W = rst->window->vertexEdges[indices[c++]];
		for (int i = 0; i < nrE; i++) E.push_back(rst->model->edgeVector[indices[c++]]);
		for (int i = 0; i < nrVe; i++) Ve.push_back(rst->model->vertices[indices[c++]]);
		//if (nrVt) Vt = rst->model->vertices[indices[i++]];
		//if (nrT) T = rst->model->edgeVector[indices[i++]];
		constructRays();
	};

	void findPrimsForESLs() {
		for (Line4& r : rays) {
			if (checkRay(r)) {
				orderDepthIntersection(r);
				std::set<int> prims;
				if (castRay(r, prims)) {
					r.prims = prims;
					putIdsInGenerators();
				}
			}
		}
	}

	void putIdsInGenerators() {
		for (Edge* e : E) e->generator_esls.insert(id);
		for (Vertex* v : Ve) v->generator_esls.insert(id);
	}

	bool checkRay(Line4 & ray) {
		if (!checkSilhouettes(ray)) return false;
		if (!nrVw) {
			if (!checkWindow(ray)) return false;
		}
		if (!checkEdgeNotVertex(ray)) return false;
		if (nrE >= 3) {
			for (Edge* e : E) {
				if (!e->intersectsRay(ray)) return false;
			}
		}
		return true;
	}

	bool checkEdgeNotVertex(Line4& ray) {
		for (Edge* e : E) {
			for (Vertex* v : e->vertices) {
				if (ray.throughPoint(v->pos, 1E-8)) return false;
			}
		}
		return true;
	}

	bool checkWindow(Line4& ray) {
		double margin = 0;
		if (nrW) margin = 1E-6;
		return rst->window->inBounds(ray, margin);
	}

	void checkInPlane() {
		if (checkInPlaneByConstruction()) {
			for (Line4& r : rays) r.inplane = true;
		}

	}

	bool checkInPlaneByConstruction() {
		if (nrVe && nrE) {
			for (Edge* e : E) {
				for (Vertex* v : Ve) {
					if (e->isPrimitive(v->triangles)) return true;
				}
			}
		}
		else if (nrE) {
			for (int i = 0; i < E.size(); i++) {
				for (int j = i + 1; j < E.size(); j++) {
					if (E[i]->isPrimitive(E[j]->triangles)) return true;
				}
			}
		}
		else if (nrVe == 2) {

		}
		return false;
	}

	void constructRays() {
		if (nrVertices == 2) {
			if (nrVw && nrVe) rays = { createVVLine(Vw, Ve[0]) };
			else if (nrVe == 2) rays = { createVVLine(Ve[0], Ve[1]) };
			//else if (nrVe && nrVt) ray1 = createVVLine(Ve[0], Vt);
		}
		else if (nrVertices == 1) {
			//if (nrVw && nrE && nrT) ray1 = createEEVLine(E[0], T, Vw);
			if (nrVw && nrE == 2) rays = { createEEVLine(E[0], E[1], Vw) };
			//if (nrW && nrVe && nrT) ray1 = createEEVLine(W, T, Ve[0]);
			else if (nrW && nrE && nrVe) rays = { createEEVLine(W, E[0], Ve[0]) };
			//if (nrW && nrE && nrVt) ray1 = createEEVLine(W, E[0], Vt);
			else if (nrVe && nrE == 2) rays = { createEEVLine(E[0], E[1], Ve[0]) };
		}
		else {
			if (nrW && nrE == 3) rays = createEEEELine(W, E);
			else if (nrE == 4) rays = createEEEELine(E);
		}
	}

	void orderDepthIntersection(Line4& ray) {
		// order the silhouettes
		std::map<double, std::set<int>> intersections;
		for (Edge* e: E) intersections[ray.depthToIntersectionWithRay(e->ray)] = e->getPrimIds();
		for (Vertex* v: Ve) intersections[ray.depthToPointOnRay(v->pos)] = getPrimIdsForVertex(v);
		//if (nrT) intersections[ray.depthToIntersectionWithRay(T->ray)] = T->getPrimIds();
		//if (nrVt) intersections[ray.depthToPointOnRay(Vt->pos)] = getPrimIdsForVertex(Vt);
		for (auto inter : intersections) {
			depth.push_back(inter.first);
			intprims.push_back(inter.second);
		}
	}

	bool vertexSilhouetteForRay(Vertex* v, Ray& ray) {
		int count = 0;
		for (Edge* e : v->edges) {
			if (e->isSilhouetteForRay(ray)) count++;
			if (count == 2) return true;
		}
		return false;
	}


	std::set<int> getPrimIdsForVertex(Vertex* v) {
		std::set<int> prims;
		for (Primitive* p : v->triangles) prims.insert(p->id);
		return prims;
	}

	bool checkSilhouettes(Line4& ray) {
		for (Edge* e : E) if (!e->isSilhouetteForRay(ray)) return false;
		for (Vertex* v : Ve) if (!vertexSilhouetteForRay(v, ray)) return false;
		return true;
	}

	//bool checkOrder(Line4& ray) {
	//	if (nrT && ray.depthToIntersectionWithRay(T->ray) != depth[depth.size() - 1]) return false;
	//	if (nrVt && ray.depthToPointOnRay(Vt->pos) != depth[depth.size() - 1]) return false;
	//	return true;
	//}

	bool castRay(Line4& ray, std::set<int>& prims) {
		double raytraceDepth; int raytracePrim;
		double cast_offset = 0;
		bool hit = false;
		double check_offset = 1E-7;
		Ray castRay(ray);

		for (int i = 0; i < depth.size(); i++) {
			hit = rst->model->getIntersectionNoAcceleration(castRay, raytracePrim, raytraceDepth);
			if (depth[i] - cast_offset + check_offset > raytraceDepth && !intprims[i].count(raytracePrim)) return false;
			else if (depth[i] - cast_offset - check_offset > raytraceDepth) return false;
			castRay = ray;
			castRay.offsetByDepth(depth[i] + 0.001);
			cast_offset = depth[i] + 0.001;
		}
		for (int p : intprims[intprims.size()-1]) {
			if (glm::dot(rst->model->triangles[p]->normal, ray.direction) < 0) prims.insert(p);
		}
		hit = rst->model->getIntersectionNoAcceleration(castRay, raytracePrim, raytraceDepth);
		if (hit) {
			prims.insert(raytracePrim);
			rst->model->triangles[raytracePrim]->esls.insert(id);
		}
		return true;
	}

	//bool checkDepth(Line4& ray) {
	//	double offset = 1E-7;
	//	for (double d : depth) {
	//		double raytraceDepth; int raytracePrim;
	//		bool hit = rst->model->getIntersectionNoAcceleration(ray, raytracePrim, raytraceDepth);
	//		if (d - offset > raytraceDepth) return false;
	//		ray.offsetByDepth(d + 0.001);
	//		offset += (d + 0.001);
	//	}
	//	return true;
	//}


	void offsetAndDirectRay(Line4& ray) {
		if (glm::dot(ray.direction, rst->window->normal) > 0) ray.inverseDir();
		double depth = rst->window->rayIntersectionDepth(ray);
		ray.offsetByDepth(depth);
	}

	std::vector<Line4> createEEEELine(std::vector<glm::dvec3> e1, std::vector<Edge*> edges) {
		std::vector<Ray> lines;
		for (Edge* e: edges) lines.push_back(e->ray);
		lines.push_back(Ray(e1[0], e1[1]));
		std::vector<Line4> lines4 = Lines4Finder::find(lines);
		for (Line4& r : lines4) {
			r.get3DfromPlucker();
			offsetAndDirectRay(r);
		}		
		return lines4;
	}

	std::vector<Line4> createEEEELine(std::vector<Edge*> edges) {
		std::vector<Ray> lines;
		for (Edge* e : edges) lines.push_back(e->ray);
		std::vector<Line4> lines4 = Lines4Finder::find(lines);
		for (Line4& r : lines4) {
			r.get3DfromPlucker();
			offsetAndDirectRay(r);
		}
		return lines4;
	}

	Line4 createEEVLine(std::vector<glm::dvec3> e1, Edge* e2, Vertex* v) {
		TriWedge e1v(v->pos, e1[0], e1[1], glm::vec3(0, 0, 0));
		glm::dvec3 from = e1v.rayIntersection(e2->ray);
		Line4 ray(from, v->pos);
		offsetAndDirectRay(ray);
		return ray;
	}


	Line4 createEEVLine(Edge* e1, Edge* e2, glm::dvec3 v) {
		TriWedge e1v(v, e1->vertices[0]->pos, e1->vertices[1]->pos, glm::vec3(0, 0, 0));
		glm::dvec3 to = e1v.rayIntersection(e2->ray);
		return Line4(v, to);
	}

	Line4 createEEVLine(Edge* e1, Edge* e2, Vertex* v) {
		TriWedge e1v(v->pos, e1->vertices[0]->pos, e1->vertices[1]->pos, glm::vec3(0, 0, 0));
		glm::dvec3 from = v->pos;
		glm::dvec3 to = e1v.rayIntersection(e2->ray);
		Line4 ray(from, to);
		offsetAndDirectRay(ray);
		return ray;
	}

	Line4 createVVLine(Vertex* v1, Vertex* v2) {
		Line4 ray(v1->pos, v2->pos);
		offsetAndDirectRay(ray);
		return ray;
	}

	Line4 createVVLine(glm::dvec3 v1, Vertex* v2) {
		return Line4(v1, v2->pos);
	}

};