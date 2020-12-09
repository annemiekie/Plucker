#ifndef CUBE_H
#define CUBE_H

#include <GL/glew.h>

#include "ray.h"
#include <glm/vec3.hpp>
#include <glm/glm.hpp>

class Cube {
    public:
    //glm::vec3 minvec, maxvec;
    glm::vec3 bounds[2];
    float size;

    Cube() {};

    Cube(glm::vec3 minvec, float size) : size(size) {
        bounds[0] = minvec;
        bounds[1] = minvec + size;
    };

    Cube(glm::vec3 minvec, glm::vec3 maxvec) {
        bounds[0] = minvec;
        bounds[1] = maxvec;
    };

    bool intersect(const Ray& r, glm::vec3& start, glm::vec3& end, glm::vec3 maindir, glm::vec3& tm) const
    {
        bool sign[3];
        sign[0] = (r.invdir.x < 0);
        sign[1] = (r.invdir.y < 0);
        sign[2] = (r.invdir.z < 0);
        float tmin, tmax, txmin, txmax, tymin, tymax, tzmin, tzmax;

        txmin = (bounds[sign[0]].x - r.origin.x) * r.invdir.x;
        txmax = (bounds[1 - sign[0]].x - r.origin.x) * r.invdir.x;
        tmin = txmin;
        tmax = txmax;
        tymin = (bounds[sign[1]].y - r.origin.y) * r.invdir.y;
        tymax = (bounds[1 - sign[1]].y - r.origin.y) * r.invdir.y;

        if ((tmin > tymax) || (tymin > tmax))
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;

        tzmin = (bounds[sign[2]].z - r.origin.z) * r.invdir.z;
        tzmax = (bounds[1 - sign[2]].z - r.origin.z) * r.invdir.z;

        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;

        //start = r.origin + std::min(tmin, tmax) * r.direction;
        //end = r.origin + std::max(tmin,tmax) * r.direction;

        glm::vec3 tmin_xyz = glm::vec3(txmin, tymin, tzmin);
        glm::vec3 tmax_xyz = glm::vec3(txmax, tymax, tzmax);
        //glm::vec3 x;
        //if (glm::dot(maindir, glm::vec3(1)) < 0) {
        float x = glm::dot(maindir * maindir, tmin_xyz);
        start = r.origin + tmin * r.direction;
        end = r.origin + tmax * r.direction;
        //}
        //else {
        //    x = maindir * maindir * tmax_xyz;
        //    start = r.origin + tmax * r.direction;
        //    end = r.origin + tmin * r.direction;
        //}
        tm = r.origin + x * r.direction;

        return true;
    };

    GLuint cubelineGeneration() {
        std::vector<GLfloat> lineSegments = {
            bounds[0].x, bounds[0].y, bounds[0].z,
            bounds[1].x, bounds[0].y, bounds[0].z,
            bounds[0].x, bounds[0].y, bounds[0].z,
            bounds[0].x, bounds[0].y, bounds[1].z,
            bounds[1].x, bounds[0].y, bounds[0].z,
            bounds[1].x, bounds[0].y, bounds[1].z,
            bounds[0].x, bounds[0].y, bounds[1].z,
            bounds[1].x, bounds[0].y, bounds[1].z,

            bounds[0].x, bounds[1].y, bounds[0].z,
            bounds[1].x, bounds[1].y, bounds[0].z,
            bounds[0].x, bounds[1].y, bounds[0].z,
            bounds[0].x, bounds[1].y, bounds[1].z,
            bounds[1].x, bounds[1].y, bounds[0].z,
            bounds[1].x, bounds[1].y, bounds[1].z,
            bounds[0].x, bounds[1].y, bounds[1].z,
            bounds[1].x, bounds[1].y, bounds[1].z,

            bounds[0].x, bounds[0].y, bounds[0].z,
            bounds[0].x, bounds[1].y, bounds[0].z,
            bounds[0].x, bounds[0].y, bounds[1].z,
            bounds[0].x, bounds[1].y, bounds[1].z,
            bounds[1].x, bounds[0].y, bounds[1].z,
            bounds[1].x, bounds[1].y, bounds[1].z,
            bounds[1].x, bounds[0].y, bounds[0].z,
            bounds[1].x, bounds[1].y, bounds[0].z,
        };

        GLuint lineVAO, lineVBO;
        glGenVertexArrays(1, &lineVAO);
        glGenBuffers(1, &lineVBO);
        glBindVertexArray(lineVAO);
        glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * lineSegments.size(), &lineSegments[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
        return lineVAO;
    };

};

#endif