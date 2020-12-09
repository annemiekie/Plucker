#pragma once

#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <CImg.h>
#include "perscamera.h"
#include "orthocamera.h"
#include "rst.h"
#include "model.h"

namespace RayTracer {

	bool intersection(int i, std::vector<Vertex>& v, Ray& r) {
		Ray edge;
		edge = Ray(v[i + 1].pos, v[i].pos);
		if (edge.side(r)) return false;
		edge = Ray(v[i + 2].pos, v[i + 1].pos);
		if (edge.side(r)) return false;
		edge = Ray(v[i].pos, v[i + 2].pos);
		if (edge.side(r)) return false;
		return true;
	};

	void imageTracer(Camera& cam, RaySpaceTree* rst, glm::ivec2 resolution, Model& model) {
		cimg_library::CImg<unsigned char> image(resolution.x, resolution.y, 1, 3, 0);
		const unsigned char color[] = { 0 , 0, 255 };

		for (int y = 0; y < resolution.y; y++) {
			int yind = (resolution.y - 1 - y);
			for (int x = 0; x < resolution.x; x++) {
				const glm::vec2 pixpos{ (x + .5f) / resolution.x * 2.0f - 1.0f, 1.0f - (y + .5f) / resolution.y * 2.0f };
				Ray ray = cam.pixRayDirection(pixpos);
				Node* leaf = rst->descend(ray);
				for (int index : leaf->primitiveSet) {
				//for (int index=0; index<vertices.size(); index+=3) {
					//if (intersection(index, vertices, ray)) {
					if (intersection(index*3, model.vertices, ray)) {
						image.draw_point(x, y, color);
						break;
					}
				}

			}
		}
		image.save("raytracer.bmp");
	};


}