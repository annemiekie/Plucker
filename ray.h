#pragma once
#define GLEW_STATIC
#include <GL/glew.h>

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <Eigen/dense>
#include <iostream>

struct Ray {
    glm::vec3 origin{ 0.0f };
    glm::vec3 direction{ 0.0f, 0.0f, -1.0f };
    glm::vec3 invdir;
    glm::vec3 u;
    glm::vec3 v;

    Ray() {};

    ~Ray() {};

    Ray(glm::vec3 q, glm::vec3 p) : origin(p) {
        u = q - p;
        v = glm::cross(q, p);
        direction = glm::normalize(u);
        invdir = 1.f / direction;
        normalize();
    };

    Ray(Eigen::VectorXd plucker, bool reverse) {
        glm::vec3 secondhalf = glm::vec3(plucker[3], plucker[4], plucker[5]);
        glm::vec3 firsthalf = glm::vec3(plucker[0], plucker[1], plucker[2]);
        u = reverse ? secondhalf : firsthalf;
        v = reverse ? firsthalf : secondhalf;
    };

    void get3DfromPlucker() {
        //glm::vec3 plane1a = glm::vec3(1, 0, 0);
        //float plane1a0 = -3.f;
        //glm::vec3 plane2a = glm::vec3(1, 1, 0);
        //float plane2a0 = 2.f;
        //glm::vec3 q = (glm::cross(plane1a, v) - plane1a0 * u) / glm::dot(plane1a, u);
        //glm::vec3 p = (glm::cross(plane2a, v) - plane2a0 * u) / glm::dot(plane2a, u);

        glm::vec3 ucrossv = glm::cross(v, u);
        float udotu = glm::dot(u, u);
        glm::vec3 p = ucrossv / udotu;
       // glm::vec3 q = p + u;
       // p -= (2.f * u);
       // q += (2.f * u);

        origin = p;
        direction = glm::normalize(u);
        invdir = 1.f / direction;

    }

    void dividefirst() {
        v /= u.x;// glm::dot(u, u);
        u /= u.x;// glm::dot(u, u);
    }

    void normalize() {
        float norm = sqrt(glm::dot(v, v) + glm::dot(u, u));
        u = u / norm;
        v = v / norm;
    }

    std::vector<double> plucker() {
        std::vector<double> pluck = { u[0], u[1], u[2], v[0], v[1], v[2] };
        return pluck;
    }

    bool side(Ray& oRay) {
        return glm::dot(u, oRay.v) + glm::dot(v, oRay.u) < 0;
    };

    float sideVal(Ray& oRay) {
        return float(glm::dot(u, oRay.v) + glm::dot(v, oRay.u));
    };

    bool sideD(Ray& oRay) {
        return glm::dot((glm::dvec3)u, (glm::dvec3)oRay.v) + glm::dot((glm::dvec3)v, (glm::dvec3)oRay.u) < 0;
    };

    float sideValD(Ray& oRay) {
        return float(glm::dot((glm::dvec3)u, (glm::dvec3)oRay.v) + glm::dot((glm::dvec3)v, (glm::dvec3)oRay.u));
    };

    void print() {
        std::cout << "U: (" << u.x << "," << u.y << "," << u.z << ") V: (" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
    }

    bool equal(Ray o, float eps) {
        dividefirst();
        o.dividefirst();
        glm::vec3 diffu = glm::abs(u - o.u);
        glm::vec3 diffv = glm::abs(v - o.v);
        if (!glm::all(glm::lessThan(diffu, glm::vec3(eps)))) return false;
        if (!glm::all(glm::lessThan(diffv, glm::vec3(eps)))) return false;
        return true;
    }
};