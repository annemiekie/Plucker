#pragma once
#define GLEW_STATIC
#include <GL/glew.h>

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>

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
    };

    bool side(Ray& oRay) {
        return glm::dot(u, oRay.v) + glm::dot(v, oRay.u) < 0;
    };

    float sideVal(Ray& oRay) {
        return float(glm::dot(u, oRay.v) + glm::dot(v, oRay.u));
    };


};