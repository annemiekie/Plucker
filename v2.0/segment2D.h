#pragma once
#include <glm/glm.hpp>
#include <glm/common.hpp>

struct Segment2D {
	glm::dvec2 p1;
	glm::dvec2 p2;

    bool intersectSegmSegm(Segment2D& o) {
        //float d = (p2.x - p1.x) * (o.p2.y - o.p1.y) - (p2.y - p1.y) * (o.p2.x - o.p1.x);
        //if (d == 0) return false;
        //float q = (p1.y - o.p1.y) * (o.p2.x - o.p1.x) - (p1.x - o.p1.x) * (o.p2.y - o.p1.y);
        //float r = q / d;
        //q = (p1.y - o.p1.y) * (o.p2.x - p1.x) - (p1.x - o.p1.x) * (p2.y - p1.y);
        //float s = q / d;
        //if (r < 0 || r > 1 || s < 0 || s > 1) return false;
        //return true;
        double s1_x = p2.x - p1.x;     double s1_y = p2.y - p1.y;
        double s2_x = o.p2.x - o.p1.x;     double s2_y = o.p2.y - o.p1.y;

        double s, t;
        s = (-s1_y * (p1.x - o.p1.x) + s1_x * (p1.y - o.p1.y)) / (-s2_x * s1_y + s1_x * s2_y);
        t = (s2_x * (p1.y - o.p1.y) - s2_y * (p1.x - o.p1.x)) / (-s2_x * s1_y + s1_x * s2_y);

        return (s >= 0 && s <= 1 && t >= 0 && t <= 1);
    }

};