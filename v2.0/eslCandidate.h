#pragma once
#include "lineThroughFourLines.h"
#include "edge.h"
#include "vertex.h"

struct ESLCandidate {
	Line4 ray;
	std::vector<Ray> lines4;
	std::vector<Split> splittingLines;
	std::vector<Vertex*> silhouetteVertices;
	std::vector<Edge*> silhouetteEdges;
	std::vector<Edge*> triangleEdges;
	std::vector<uint64_t> indices;
	bool inBox = true;
	
	void fillUpLines() {
		for (Split& s : splittingLines) lines4.push_back(s.ray);
		for (Vertex* v : silhouetteVertices) {
			lines4.push_back(v->edges[0]->ray);
			lines4.push_back(v->edges[1]->ray);
		}
		for (Edge* e : silhouetteEdges) lines4.push_back(e->ray);
		for (Edge* e : triangleEdges) lines4.push_back(e->ray);
	}
};