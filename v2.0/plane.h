#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include "ray.h"
#include <vector>
class Plane {
public:
    glm::dvec3 normal;
    double constant = 0;

    Plane() {};

    Plane(Plane* plane) {
        normal = plane->normal;
        constant = plane->constant;
    };

    Plane(glm::dvec3 n, double c) : normal(n), constant(c) {};

    Plane(glm::dvec3 point, glm::dvec3 normal) : normal(normal) {
        constant = glm::dot(normal, point);
    }

    Plane(std::vector<glm::dvec3> vertices, glm::dvec3 voppo) {
        glm::dvec3 e1 = vertices[1] - vertices[0];
        glm::dvec3 e2 = vertices[2] - vertices[0];
        normal = glm::normalize(glm::cross(e1, e2));
        constant = glm::dot(normal, vertices[0]);
        if (glm::dot(normal, voppo) - constant > 0) {
            normal *= -1;
            constant *= -1;
        }
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

    Plane(Ray &r1, Ray &r2) {
        normal = glm::normalize(glm::cross(r1.direction, r2.direction));
        constant = glm::dot(normal, r1.origin);
    }

    Plane(glm::dvec3& dir1, glm::dvec3& dir2, glm::dvec3 point) {
        normal = glm::normalize(glm::cross(dir1, dir2));
        constant = glm::dot(normal, point);
    }

    bool rayInPlane(const Ray& r, double eps = 1E-10) {
        return (pointOnPlane(r.origin, eps) && pointOnPlane(r.origin + r.direction, eps)); //glm::dot(r.direction, normal) < eps && pointOnPlane(r.origin, eps);
    }

    double distToPoint(glm::dvec3 pt) {
        return glm::dot(normal, pt) - constant;
    }

    bool pointOnPlane(glm::dvec3 pt, double eps = 1E-10) {
        return  fabs(distToPoint(pt)) < eps;
    }

    bool pointOnPositiveSide(glm::dvec3 pt) {
        return distToPoint(pt) > 0;
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

    bool equal(Plane* other, double thres = 1E-10) {
        if (glm::length(glm::normalize(normal) - glm::normalize(other->normal)) > thres) return false;
        if (fabs(constant - other->constant) > thres) return false;
        return true;
    }

    bool isInFrontOf(glm::dvec3 v1, glm::dvec3 v2) {
        if (glm::length(v1-v2)<1E-10) return false;
        return fabs(distToPoint(v1)) < fabs(distToPoint(v2));
    }

    bool isInFrontOf(glm::dvec3 v, std::vector<glm::dvec3>& e) {
        return  isInFrontOf(v, e[0]) || isInFrontOf(v, e[1]);
    }

    bool isInFrontOf(std::vector<glm::dvec3>& e, glm::dvec3 v) {
        return isInFrontOf(e[0], v) || isInFrontOf(e[1], v);
    }

    bool isInFrontOf(std::vector<glm::dvec3>& e1, std::vector<glm::dvec3>& e2)  {
        return  isInFrontOf(e1, e2[0]) || isInFrontOf(e1, e2[1]);
    }


    ~Plane() {};

};