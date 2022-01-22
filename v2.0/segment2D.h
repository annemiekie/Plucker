#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>

struct Segment2D {
	glm::vec2 p1;
	glm::vec2 p2;

    bool intersectSegmSegm(Segment2D o) {
        float d = (p2.x - p1.x) * (o.p2.y - o.p1.y) - (p2.y - p1.y) * (o.p2.x - o.p1.x);
        if (d == 0) return false;
        float q = (p1.y - o.p1.y) * (o.p2.x - o.p1.x) - (p1.x - o.p1.x) * (o.p2.y - o.p1.y);
        float r = q / d;
        q = (p1.y - o.p1.y) * (o.p2.x - p1.x) - (p1.x - o.p1.x) * (p2.y - p1.y);
        float s = q / d;
        if (r < 0 || r > 1 || s < 0 || s > 1) return false;
        return true;
    }

};