#ifndef EDGE_H
#define EDGE_H

#pragma once
struct Edge {
	std::set<int> vertices;
	//int v1;
	//int v2;
	mutable std::vector<int> triangles = std::vector<int>();
	mutable int index;
};

struct cmp_by_v {
	bool operator()(const Edge& a, const Edge& b) const {
		return a.vertices < b.vertices;
		//return a.v1 < b.v1;
	}
};

#endif