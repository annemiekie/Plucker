#pragma once

#include "ray.h"
#include "sample.h"
#include <utility> 
#include <sstream>
#include <set>
#include <vector>

class Node {
	public:
	bool leaf = false;
	Node *leftNode = nullptr;
	Node *rightNode = nullptr;
	Node* parent = nullptr;
	Ray splitter = Ray();
	std::set<int> primitiveSet = std::set<int>();
	std::vector<Sample> primAndRayVector = std::vector<Sample>();
	int index = 0;
	int depth = 0;

	Node() : leaf(false) {};

	~Node() {};

	Node(int ind, int depth) : index(ind), depth(depth), leaf(false) { };

	void insert(int pri) {
		primitiveSet.insert(pri);
	}

	void insert(int pri, Ray ray) {
		primitiveSet.insert(pri);
		primAndRayVector.push_back(Sample{ray, pri});
	};
};
