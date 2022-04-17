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
	std::vector<Edge*> allEdges;
	bool inBox = true;
	
	void fillUpLines() {
		for (SplitSide& s : splittingLines) {
			lines4.push_back(s.ray);
			if (s.edge != NULL) allEdges.push_back(s.edge);
		}
		for (Vertex* v : silhouetteVertices) {
			lines4.push_back(v->edges[0]->ray);
			lines4.push_back(v->edges[1]->ray);
			for (Edge* ve : v->edges) allEdges.push_back(ve);
		}
		for (Edge* e : silhouetteEdges) {
			lines4.push_back(e->ray);
			allEdges.push_back(e);
		}
		for (Edge* e : triangleEdges) {
			lines4.push_back(e->ray);
			allEdges.push_back(e);
		}
	}

	bool combinationDegenerate() {
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

		for (Edge* se : silhouetteEdges) {
			for (Vertex* v1 : se->vertices) {
				for (Edge* te : triangleEdges) {
					for (Vertex* v2 : te->vertices) {
						if (v1->id == v2->id) return true;
					}
				}
			}
		}

		for (int x = 0; x < lines4.size(); x++) {
			for (int y = x + 1; y < lines4.size(); y++) {
				if (lines4[x].equal(lines4[y])) return true;
			}
		}

		for (Vertex* v : triangleVertex) {
			for (Edge* e : silhouetteEdges) {
				if (e->ray.throughVertex(v)) return true;
			}
			for (SplitSide& s : splittingLines) {
				if (s.edge != NULL) {
					if (s.edge->ray.throughVertex(v)) return true;
				}
			}
		}
		for (Vertex* v : silhouetteVertices) {
			for (Edge* e : triangleEdges) {
				if (e->ray.throughVertex(v)) return true;
			}
			for (SplitSide& s : splittingLines) {
				if (s.edge != NULL) {
					if (s.edge->ray.throughVertex(v)) return true;
				}
			}
		}
		return false;
	}
};