#pragma once
#define GLEW_STATIC
#include <GL/glew.h>

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <Eigen/dense>
#include <iostream>

struct Ray {
    glm::dvec3 origin{ 0.0f };
    glm::dvec3 direction{ 0.0f, 0.0f, -1.0f };
    glm::dvec3 invdir;
    glm::dvec3 u;
    glm::dvec3 v;

    Ray() {};

    ~Ray() {};

    Ray(glm::dvec3 q, glm::dvec3 p) : origin(p) {
        u = q - p;
        v = glm::cross(q, p);
        direction = glm::normalize(u);
        invdir = 1. / direction;
        normalize();
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

    //void normalizeU() {
    //    double norm = sqrt(glm::dot(u, u));
    //    u = u / norm;
    //    v = v / norm;
    //}
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

    void print() {
        std::cout << "U: (" << u.x << "," << u.y << "," << u.z << ") V: (" << v.x << "," << v.y << "," << v.z << ")";
    }

    bool equal(Ray o, double eps) {
        normalize();
        o.normalize();
        glm::dvec3 diffv = glm::abs(v/u.x - o.v/o.u.x);
        glm::dvec3 diffu = glm::abs(u/u.x - o.u/o.u.x);
        if (!glm::all(glm::lessThan(diffu, glm::dvec3(eps)))) return false;
        if (!glm::all(glm::lessThan(diffv, glm::dvec3(eps)))) return false;
        return true;
    }
};