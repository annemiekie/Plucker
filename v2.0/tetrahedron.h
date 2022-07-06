#pragma once
#include "triWedge.h"
#include <vector>

class Tetrahedron {
public:
	std::vector<TriWedge> wedges;
	bool unbound = false;
	std::vector<std::vector<double>> m = { {1,1,1,1}, {-1,-1, 1,1}, {1,1,-1,-1} };

	Tetrahedron() {};
	~Tetrahedron() {};

	Tetrahedron(std::vector<TriWedge> wedges, bool unbound = true) : wedges(wedges), unbound(unbound) {};

	Tetrahedron(std::vector<glm::dvec3>& e1v, std::vector<glm::dvec3>& e2v, bool unbound = true) : unbound(unbound) {
		wedges = {
			TriWedge(e2v[0], e1v[0], e1v[1], e2v[1]),
			TriWedge(e2v[1], e1v[0], e1v[1], e2v[0]),
			TriWedge(e1v[0], e2v[0], e2v[1], e1v[1]),
			TriWedge(e1v[1], e2v[0], e2v[1], e1v[0]),
		};
	}

	bool segmInTetra(glm::dvec3 v1, glm::dvec3 v2, double err = 0) {
		if (pointInTetra(v1, 0) || pointInTetra(v2, 0)) return true;
		glm::dvec3 intPt;
		for (TriWedge& w : wedges) if (w.raySegmentIntersection(v1, v2, intPt)) return true;
			//if (pointInTetra(intPt, 1E-8)) return true;
		return false;
	}

	bool pointInTetra(glm::dvec3 pt, double err) {
		int volNr = unbound ? 3 : 1;
		for (int i=0; i<volNr; i++) {
			bool inVolume = true;
			for (int j = 0; j < wedges.size(); j++) {
				double x = glm::dot(wedges[j].normal * m[i][j], pt) - wedges[j].constant * m[i][j];
				if (glm::dot(wedges[j].normal * m[i][j], pt) - wedges[j].constant * m[i][j] > err) {
					inVolume = false;
					break;
				}
			}
			if (inVolume) return true;
		}
		return false;
	}

};