#pragma once
#define GLEW_STATIC
#include <GL/glew.h>

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <Eigen/dense>
#include <iostream>
#include "vertex.h"

class Ray {
public:
    glm::dvec3 origin{ 0.0f };
    glm::dvec3 direction{ 0.0f, 0.0f, -1.0f };
    glm::dvec3 invdir;
    glm::dvec3 u;
    glm::dvec3 v;
    int index = -1;
   // bool checkboth = true;

    Ray() {};

    ~Ray() {};

    Ray(const Ray& r) {
        origin = r.origin;
        direction = r.direction;
        invdir = r.invdir;
        u = r.u;
        v = r.v;
        index = r.index;
    }

    Ray(glm::dvec3 from, glm::dvec3 to, int index = -1) : origin(from), index(index) {
        u = to - from;
        v = glm::cross(to, from);
        direction = glm::normalize(u);
        invdir = 1. / direction;
        normalize();;
    };

    Ray(Eigen::VectorXd plucker, bool reverse) {
        glm::dvec3 secondhalf = glm::dvec3(plucker[3], plucker[4], plucker[5]);
        glm::dvec3 firsthalf = glm::dvec3(plucker[0], plucker[1], plucker[2]);
        u = reverse ? secondhalf : firsthalf;
        v = reverse ? firsthalf : secondhalf;
    };

    void get3DfromPlucker() {
        glm::dvec3 ucrossv = glm::cross(v, u);
        double udotu = glm::dot(u, u);
        glm::dvec3 p = ucrossv / udotu;
        origin = p;
        direction = glm::normalize(u);
        invdir = 1. / direction;
    }

    void dividefirst() {
        v /= u.x;// glm::dot(u, u);
        u /= u.x;// glm::dot(u, u);
    }

    void inverseDir() {
        direction = -direction;
        invdir = 1. / direction;
        u = -u;
        v = -v;
    }

    void normalize() {
        double norm = sqrt(glm::dot(v, v) + glm::dot(u, u));
        u = u / norm;
        v = v / norm;
    }

    std::vector<double> plucker() {
        std::vector<double> pluck = { u[0], u[1], u[2], v[0], v[1], v[2] };
        return pluck;
    }

    bool side(const Ray& oRay) const {
        return glm::dot(u, oRay.v) + glm::dot(v, oRay.u) < 0;
    };

    double sideVal(const Ray& oRay) const {
        return glm::dot(u, oRay.v) + glm::dot(v, oRay.u);
    };

    bool intersect(const Ray& oRay, double thres = 1E-8) const {
        if (fabs(sideVal(oRay)) < thres) return true;
        return false;
    };

    bool intersectsWithRayAtDepth(Ray& oRay, double depth, float thres = 1E-8) {
        glm::dvec3 intersectionTs = (origin + direction * depth - oRay.origin) / oRay.direction;
        return (fabsf(intersectionTs.x - intersectionTs.y) < thres && fabsf(intersectionTs.x - intersectionTs.z) < thres);
    }

    glm::dvec3 pointOfintersectWithRay(const Ray& oRay) const {
        glm::dvec3 N(0, 0, 1);
        return (glm::dot(v, N) * oRay.u - glm::dot(oRay.v, N) * u - glm::dot(v, oRay.u) * N) /
                   glm::dot(glm::cross(u, oRay.u), N);
    }

    void print() {
        std::cout << "U: (" << u.x << "," << u.y << "," << u.z << ") V: (" << v.x << "," << v.y << "," << v.z << ")";
    }

    void printod() {
        std::cout << "Origin: (" << origin.x << "," << origin.y << "," << origin.z << ") Direction: (" << direction.x << "," << direction.y << "," << direction.z << ")";
    }

    bool throughVertex(Vertex* vertex, double eps = 1E-8) {
        return throughPoint(vertex->pos, eps);
        //glm::dvec3 dir = vertex - origin;
        //if (dir.x * direction.x < 0) dir = -dir;
        //double dist = glm::distance(glm::normalize(dir), glm::normalize(direction));
        //return (dist < eps);
    }

    void offsetByDepth(double depth) {
        origin += (depth * direction);
    }

    double pointToRayDist(glm::dvec3 pt) {
        return glm::length(glm::cross(pt - origin, pt - origin + direction)) /
                        glm::length(direction);

    }

    bool throughPoint(glm::dvec3 pt, double eps = 1E-8) {
        return (pointToRayDist(pt) < eps);
    }

    double depthToIntersectionWithRay(Ray& oRay) const {
        glm::dvec3 intersect = pointOfintersectWithRay(oRay);
        return depthToPointOnRay(intersect);
    }

    double depthToPointOnRay(glm::dvec3 pt) const {
        return (pt.x - origin.x) / (direction.x);
    }

    //bool inPlane(Ray& ray1, Ray& ray2, double eps = 1E-8) {
    //    glm::dvec3 cross = glm::cross(glm::normalize(ray1.direction), glm::normalize(ray2.direction));
    //    double v = fabsf(glm::dot(direction, cross));
    //    return (v < eps);
    //}

    // Checks whether rays are the same
    // Not checking orientation 
    bool equal(Ray o, double eps = 1E-8) {
        normalize();
        o.normalize();
        if ((glm::distance(u, o.u) < eps && glm::distance(v, o.v) < eps)) return true;
        o.inverseDir();
        if ((glm::distance(u, o.u) < eps && glm::distance(v, o.v) < eps)) return true;
        return false;
        // glm::dvec3 diffv = glm::abs(v - o.v);
        // glm::dvec3 diffu = glm::abs(u - o.u);
        // if (!glm::all(glm::lessThan(diffu, glm::dvec3(eps)))) return false;
        // if (!glm::all(glm::lessThan(diffv, glm::dvec3(eps)))) return false;
        // return true;
    }
};