#pragma once

#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <CImg.h>
#include "perscamera.h"
#include "orthocamera.h"
#include "rst.h"
#include "model.h"

namespace RayTracer {

	float intersection(int i, std::vector<Vertex>& v, Ray& r, float &diff) {
		Ray edge;
		edge = Ray(v[i + 1].pos, v[i].pos);
		if (edge.side(r)) return 0;
		edge = Ray(v[i + 2].pos, v[i + 1].pos);
		if (edge.side(r)) return 0;
		edge = Ray(v[i].pos, v[i + 2].pos);
		if (edge.side(r)) return 0;

		glm::dvec3 v0 = v[i].pos;
		glm::dvec3 v1 = v[i + 1].pos;
		glm::dvec3 v2 = v[i + 2].pos;

		glm::dvec3 N = glm::cross(v2 - v0, v1 - v0);
		glm::dvec3 L = glm::vec3(-3, 3, -3);
		diff = std::max(glm::dot(-glm::normalize(N), glm::normalize(L)), 0.);
		float ndotdir = glm::dot(N, r.direction);
		float t = glm::dot(N, v0 - r.origin) / ndotdir;
		return t;
	};

	bool intersectionVal(int i, std::vector<Vertex>& v, Ray& r) {
		Ray edge;
		float edgeVal;
		edge = Ray(v[i + 1].pos, v[i].pos);
		edgeVal = edge.sideVal(r);
		std::cout << edgeVal << std::endl;
		//if (abs(edgeVal) < 1E-10) return true;
		edge = Ray(v[i + 2].pos, v[i + 1].pos);
		edgeVal = edge.sideVal(r);
		std::cout << edgeVal << std::endl;

		//if (abs(edgeVal) < 1E-10) return true;
		edge = Ray(v[i].pos, v[i + 2].pos);
		edgeVal = edge.sideVal(r);
		std::cout << edgeVal << std::endl;

		//if (abs(edgeVal) < 1E-10) return true;
		return false;
	};

	//bool intersection2(int i, std::vector<Vertex>& v, Ray& r) {
	//	glm::vec3 v0 = v[i].pos;
	//	glm::vec3 v1 = v[i+1].pos;
	//	glm::vec3 v2 = v[i+2].pos;
	//
	//	glm::vec3 N = glm::cross(v1 - v0, v2 - v0);
	//	float D = glm::dot(N, v0);
	//
	//	float ndotdir = glm::dot(N, r.direction);
	//	//if (fabs(ndotdir) < 1E-10) // almost 0 
	//	//	return false;
	//
	//	float t = -(glm::dot(N, r.origin) + D) / ndotdir;
	//	if (t < 0) return false; // the triangle is behind 
	//
	//	glm::vec3 P = r.origin + t * r.direction;
	//
	//	glm::vec3 C; // vector perpendicular to triangle's plane 
	//
	//	// edge 0
	//	glm::vec3 edge0 = v1 - v0;
	//	glm::vec3 vp0 = P - v0;
	//	C = glm::cross(edge0,vp0);
	//	if (glm::dot(N,C) < 0) return false; // P is on the right side 
	//
	//	// edge 1
	//	glm::vec3 edge1 = v2 - v1;
	//	glm::vec3 vp1 = P - v1;
	//	C = glm::cross(edge1, vp1);
	//	if (glm::dot(N, C) < 0)  return false; // P is on the right side 
	//
	//	// edge 2
	//	glm::vec3 edge2 = v0 - v2;
	//	glm::vec3 vp2 = P - v2;
	//	C = glm::cross(edge2, vp2);
	//	if (glm::dot(N, C) < 0) return false; // P is on the right side; 
	//
	//	return true;
	//}


	std::vector<int> imageTracer(Camera& cam, RaySpaceTree* rst, glm::ivec2 resolution, Model& model) {
		cimg_library::CImg<unsigned char> image(resolution.x, resolution.y, 1, 3, 0);
		std::vector<int> triNums = std::vector<int>(resolution.x*resolution.y);
		//#pragma omp parallel for
		for (int y = 0; y < resolution.y; y++) {
			std::cout << y << " ";
			for (int x = 0; x < resolution.x; x++) {
				const glm::vec2 pixpos{ (x + .5f) / resolution.x * 2.0f - 1.0f, 1.0f - (y + .5f) / resolution.y * 2.0f };
				Ray ray = cam.pixRayDirection(pixpos);
				Node* leaf = rst->descend(ray);
				float depth = 10000;
				float diff = 0.f;
				const unsigned char colorb[] = { 255,255, 255 };
				image.draw_point(x, y, colorb);
				for (int index : leaf->primitiveSet) {
				//for (int index =0; index < model.primsize; index++) {
					float t = intersection(index * 3, model.vertices, ray, diff);
					if (t >0 && t < depth) {
						char col = char(diff * 255.f);
						
						int leafnum = rst->numOfLeaf(leaf->index);
						char b = char(255.f / rst->noLeaves * leafnum);

						const unsigned char color[] = { b , 0, 255-b };
						image.draw_point(x, y, color);
						//image.draw_point(x, y, colorb);
						triNums[y * resolution.x + x] = index+1;
						depth = t;
					}
				}
			}
		}
		std::cout << std::endl;
		image.save("raytracer.bmp");
		return triNums;
	};


}