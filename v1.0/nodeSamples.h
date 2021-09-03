#pragma once
#include <vector>
#include <set>
#include "node.h"

struct nodeSamples {
	std::set<int > triangles;
	std::vector < std::pair<int, int>> samples;
	Node* node;
	int level;
};