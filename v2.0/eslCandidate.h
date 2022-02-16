#pragma once
#include "lineThroughFourLines.h"
#include "edge.h"
#include "vertex.h"

struct ESLCandidate {
	Line4 ray;
	std::vector<Ray> lines4;
	std::vector<Vertex*> silhVertices;
	std::vector<Edge*> silhEdges;
	bool inBox = false;
};