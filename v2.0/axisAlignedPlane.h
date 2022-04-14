#pragma once
#include "plane.h"

class AxisAlignedPlane : public virtual Plane {
public:
	AxisAlignedPlane() {};
	glm::dvec3 aaPosition;
	int constantCoord = 0;

	AxisAlignedPlane(AxisAlignedPlane* plane) : Plane(plane) {
		aaPosition = plane->aaPosition;
		constantCoord = plane->constantCoord;
	};

	AxisAlignedPlane(glm::dvec3 point, glm::dvec3 normal) : Plane(point, normal) {
		int count = 0;
		for (int i = 0; i < 3; i++) {
			if (normal[i] == 0) count++;
			else {
				aaPosition[i] = point[i];
				constantCoord = i;
			}
		}
		if (count != 2) std::cout << "Not axis aligned!" << std::endl;
	};

	~AxisAlignedPlane() {};


};