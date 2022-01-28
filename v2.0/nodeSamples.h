#pragma once
#include <vector>
#include <set>
#include "node.h"
#include "sample.h"
#include "edge.h"
struct NodeSamples {
	int level;
	std::vector<SampleInd> samples;
	std::set<int> triangles;
	std::set<Edge*, Edge::cmp_ptr> edges;
	Node* node;

};