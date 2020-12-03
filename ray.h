#pragma once
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

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
        float x = glm::dot(u, oRay.v) + glm::dot(v, oRay.u);
        return x < 0;
    };

    float sideVal(Ray& oRay) {
        return float(glm::dot(u, oRay.v) + glm::dot(v, oRay.u));
    };

    static Ray getRandom() {
        int t = (int)time(NULL);
        srand(t);
        float r1 = static_cast <float> (RAND_MAX) / 10.f;
        float x = float(rand()) / r1;
        float y = float(rand()) / r1;
        float z = float(rand()) / r1;
        float dx = float(rand()) / r1;
        float dy = float(rand()) / r1;
        float dz = float(rand()) / r1;
        return Ray(glm::vec3(x, y, z), glm::vec3(dx, dy, dz));
    }

};