#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include "ray.h"
class Plane {
public: 
	glm::dvec3 normal;
	float constant = 0;

    Plane() {};

    ~Plane() {};

    Plane(glm::dvec3 n, float c) : normal(n), constant(c) {};

    Plane(glm::dvec3 point, glm::dvec3 normal) : normal(normal) {
        constant = glm::dot(normal, point);
    }

    Plane(glm::dvec3 v1main, glm::dvec3 v2, glm::dvec3 v3, glm::dvec3 voppo) {
        glm::dvec3 e1 = v2 - v1main;
        glm::dvec3 e2 = v3 - v1main;
        normal = glm::normalize(glm::cross(e1, e2));
        constant = glm::dot(normal, v1main);
        if (glm::dot(normal, voppo) - constant > 0) {
            normal *= -1;
            constant *= -1;
        }
    }

    bool rayInPlane(Ray& r, float eps = 1E-10) {
        return (pointOnPlane(r.origin, eps) && pointOnPlane(r.origin + r.direction, eps));
    }

    bool pointOnPlane(glm::dvec3 pt, float eps = 1E-10) {
        return  fabsf(glm::dot(normal, pt) - constant) < eps;
    }

    bool pointOnPositiveSide(glm::dvec3 pt) {
        return (glm::dot(normal, pt) - constant) > 0;
    }

    double rayIntersectionDepth(const Ray& r) {
        return (-glm::dot(r.origin, normal) + constant) / glm::dot(r.direction, normal);
    }

    virtual glm::vec3 rayIntersection(const Ray& r) {
        return r.origin + rayIntersectionDepth(r) * r.direction;
    }

    double rayIntersection(glm::dvec3& v1, glm::dvec3& v2) {
        return (-glm::dot(v1, normal) + constant) / glm::dot(v2-v1, normal);
    }

    virtual bool raySegmentIntersection(glm::dvec3 v1, glm::dvec3 v2, glm::dvec3& pt) {
        double t = rayIntersection(v1, v2);
        if (t < 0 || t > 1) return false;
        pt = v1 + t *(v2-v1);
        return true;
    }
    
};