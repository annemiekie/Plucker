#pragma once
#include "polygon.h"

class Square : virtual public Polygone {
public:


	Square(glm::dvec3 n, float c) : Polygone(n, c) { };
	Square() {};
	~Square() {};

	Square(glm::dvec3 v1main, glm::dvec3 v2, glm::dvec3 v3, glm::dvec3 voppo) : Plane(v1main, v2, v3, voppo) {
		vertices = std::vector<glm::dvec3>(4);
		vertices[0] = v1main;
		vertices[1] = v2;
		vertices[2] = v2 + v3 - v1main;
		vertices[3] = v3;
		makeDirectedEdges();
	};

}; 

