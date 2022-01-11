#ifndef EDGE_H
#define EDGE_H

#pragma once

struct Edge {
	std::vector<int> v;
	mutable std::vector<int> triangles = std::vector<int>();
	mutable int index;
	bool operator<(const Edge& e) const {
		return  e.v[0] < this->v[0] || (e.v[0] == this->v[0] && e.v[1] < this->v[1]); 
	}
};


#endif