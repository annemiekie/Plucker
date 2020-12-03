#pragma once

#include "ray.h"
#include <utility> 
#include <sstream>
#include <set>
#include <vector>

class Node {
	public:
	bool leaf = false;
	Node *leftNode = nullptr;
	Node *rightNode = nullptr;
	Ray splitter = Ray();
	std::set<int> primitiveSet = std::set<int>();
	std::vector<std::pair<int, Ray>> primAndRayVector = std::vector<std::pair<int, Ray>>();
	int index = 0;

	Node() : leaf(false) {};

	~Node() {};

	Node(int ind) : index(ind), leaf(false) { };

	void insert(int pri, Ray ray) {
		primitiveSet.insert(pri);
		primAndRayVector.push_back(std::pair<int,Ray>(pri, ray));
	};
};
