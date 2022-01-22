#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include "ray.h"
class Plane {
public: 
	glm::vec3 normal;
	float constant = 0;

    Plane() {};

    ~Plane() {};

    Plane(glm::vec3 n, float c) : normal(n), constant(c) {};

    Plane(glm::vec3 v1main, glm::vec3 v2, glm::vec3 v3, glm::vec3 voppo) {
        glm::vec3 e1 = v2 - v1main;
        glm::vec3 e2 = v3 - v1main;
        normal = glm::normalize(glm::cross(e1, e2));
        constant = glm::dot(normal, v1main);
        if (glm::dot(normal, voppo) - constant > 0) {
            normal *= -1;
            constant *= -1;
        }
    }

    virtual glm::vec3 rayIntersection(const Ray& r) {
        double t = (-glm::dot((glm::vec3)r.origin, normal) + constant) / glm::dot((glm::vec3)r.direction, normal);
        return r.origin + t * r.direction;
    }

    float rayIntersection(glm::vec3& v1, glm::vec3& v2) {
        return (-glm::dot(v1, normal) + constant) / glm::dot(v2-v1, normal);
    }

    virtual bool raySegmentIntersection(glm::vec3 v1, glm::vec3 v2, glm::vec3& pt) {
        float t = rayIntersection(v1, v2);
        if (t < 0 || t > 1) return false;
        pt = v1 + t *(v2-v1);
        return true;
    }
    
};