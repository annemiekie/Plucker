#pragma once
#include <vector>
#include <set>
#include "node.h"
#include "sample.h"
#include "edge.h"
struct nodeSamples {
	std::set<int> triangles;
	std::vector<SampleInd> samples;
	std::set<Edge*, Edge::cmp_ptr> edges;
	Node* node;
	int level;
};