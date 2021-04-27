#ifndef EDGE_H
#define EDGE_H

#pragma once
// Per-vertex data
struct Edge {
	std::set<int> vertices;
	mutable std::vector<int> triangles = std::vector<int>();
};

struct cmp_by_v {
	bool operator()(const Edge& a, const Edge& b) const {
		return a.vertices < b.vertices;
	}
};

#endif