#pragma once
#include <vector>
#include "vertex.h"
#include "edge.h"
#include "split.h"

struct FourLines {
	std::vector<Split> splittingLines;
	std::vector<Vertex*> silhouetteVertices;
	std::vector<Edge*> silhouetteEdges;
	std::vector<Edge*> trianglesEdges;
	std::vector<Ray> lines;
};