#pragma once
#include "lineThroughFourLines.h"
#include "edge.h"
#include "vertex.h"
#include "splitSide.h"

struct ESLCandidate {
	Line4 ray;
	std::vector<Ray> lines4;
	std::vector<SplitSide> splittingLines;
	std::vector<Vertex*> silhouetteVertices;
	std::vector<Edge*> silhouetteEdges;
	std::vector<Edge*> triangleEdges;
	std::vector<Vertex*> triangleVertex;
	std::vector<uint64_t> indices;
	std::set<Edge*, Edge::cmp_ptr> allEdges;
	bool inBox = true;
	bool inPlane = false;
	bool inplaneNotVertex = false;
	bool throughVertex = false;
//	bool stabBack = false;
	
	void fillUpLines() {
		for (SplitSide& s : splittingLines) {
			lines4.push_back(s.ray);
			if (s.edge != NULL) allEdges.insert(s.edge);
			indices.push_back(s.ray.index);
		}
		for (Vertex* v : silhouetteVertices) {
			lines4.push_back(v->edges[0]->ray);
			lines4.push_back(v->edges[1]->ray);
			for (Edge* ve : v->edges) allEdges.insert(ve);
			indices.push_back(v->id);
		}
		for (Edge* e : silhouetteEdges) {
			lines4.push_back(e->ray);
			allEdges.insert(e);
			indices.push_back(e->id);
		}
		for (Edge* e : triangleEdges) {
			lines4.push_back(e->ray);
			allEdges.insert(e);
			indices.push_back(e->id);
		}
		for (Vertex* v : triangleVertex) {
			for (Edge* ve : v->edges) allEdges.insert(ve);
		}
	}

	bool doublesSplit(Split& s) {
		if (s.edge != NULL) {
			for (Edge* ae : allEdges) {
				if (s.edge->id == ae->id) return true;
				if (s.edge->ray.equal(ae->ray)) return true;
			}
		}
		else {
			for (SplitSide& sa : splittingLines) {
				if (sa.id == s.id) return true;
			}
		}
		return false;
	}

	bool doublesEdge(Edge* e) {
		for (Edge* ae : allEdges) {
			if (e->id == ae->id) return true;
			if (e->ray.equal(ae->ray)) return true;
		}
		return false;
	}

	void putRay(Ray& r) {
		ray = r;
		if (triangleVertex.size() > 0) {
			if (!ray.throughPoint(triangleVertex[0]->pos)) {
				inPlane = true;
				inPlaneNotVertex();
			}
			else throughVertex = true;
		}
	}

	void inPlaneNotVertex() {
		inplaneNotVertex = true;
		if (triangleVertex.size() == 0) return;
		triangleVertex = std::vector<Vertex*>();
		allEdges = std::set<Edge*, Edge::cmp_ptr>();
		for (Edge* e : triangleEdges) allEdges.insert(e);
		for (SplitSide& s : splittingLines) {
			if (s.edge != NULL) allEdges.insert(s.edge);
		}
		for (Vertex* v : silhouetteVertices) {
			for (Edge* ve : v->edges) allEdges.insert(ve);
		}
		for (Edge* e : silhouetteEdges) allEdges.insert(e);
	}

	bool combinationDegenerate(Primitive* prim) {
		//for (int i = 0; i < allEdges.size(); i++) {
		//	for (Vertex* v1 : allEdges[i]->vertices) {
		//		for (int j = i+1; j < allEdges.size(); j++) {
		//			for (Vertex* v2 : allEdges[j]->vertices) {
		//				if (v1->id == v2->id) return true;
		//			}
		//		}
		//		for (Vertex* v : triangleVertex) {
		//			if (v1->id == v->id) return true;
		//		}
		//	}
		//}

		for (int x = 0; x < lines4.size(); x++) {
			for (int y = x + 1; y < lines4.size(); y++) {
				if (lines4[x].equal(lines4[y])) return true;
			}
		}

		for (Vertex* v : triangleVertex) {
			for (Edge* e : silhouetteEdges) {
				if (e->ray.throughPoint(v->pos)) return true;
			}
			for (SplitSide& s : splittingLines) {
				if (s.edge != NULL) {
					if (s.edge->ray.throughPoint(v->pos)) return true;
				}
			}
		}
		for (Vertex* v : silhouetteVertices) {
			for (Edge* e : triangleEdges) {
				if (e->ray.throughPoint(v->pos)) return true;
			}
			for (SplitSide& s : splittingLines) {
				if (s.edge != NULL) {
					if (s.edge->ray.throughPoint(v->pos)) return true;
				}
			}
			if (prim->getPlane().pointOnPlane(v->pos, 1E-8)) return true;
		}

		for (Edge* e : silhouetteEdges) {
			if (prim->getPlane().rayInPlane(e->ray, 1E-8)) return true;
		}

		for (SplitSide& s : splittingLines) {
			if (s.edgeIsTriangle) return true;
			if (prim->getPlane().rayInPlane(s.ray, 1E-8)) return true;
		}


		return false;
	}

	void print() {
		for (SplitSide& s : splittingLines) {
			std::cout << "S:" << s.id;
			if (s.edge != NULL) std::cout << "(e):" << s.edge->id;
			else std::cout << "\t";
			std::cout << "\t";
		}
		for (Vertex* v : silhouetteVertices) {
			std::cout << "V:" << v->id << " (e):";
			for (Edge* ve : v->edges) std::cout << ve->id << " ";
			std::cout << "\t";
		}
		for (Edge* e : silhouetteEdges) {
			std::cout << "E:" << e->id << "\t";
		}
		if (triangleVertex.size() != 0) {
			for (Vertex* v : triangleVertex) {
				std::cout << "V:" << v->id << "(t):";
				for (Edge* ve : v->edges) std::cout << ve->id << " ";
				std::cout << "\t";
			}
		}
		else {
			for (Edge* e : triangleEdges) {
				std::cout << "T:" << e->id << " ";
				std::cout << "\t";
				if (e->id < 10) std::cout << "\t";
			}
		}
		std::cout << std::endl;
	}
};