#ifndef SPHERE_H
#define SPHERE_H

#include <GL/glew.h>

#include "ray.h"
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>
#define _USE_MATH_DEFINES
#include "math.h"

#include <vector>
#include "vertex.h"
#include "geoObject.h"

class Sphere : public GeoObject {
public:
    float radius = 0.f;
    int indSize = 0;

    Sphere() {};

    Sphere(glm::vec3 center, float radius) : GeoObject(center), radius(radius) {   };

    virtual bool intersect(const Ray &line, glm::vec3 &start, glm::vec3 &end, bool getcolor = false, glm::vec3& maindir = glm::vec3(0), glm::vec3& color = glm::vec3()) override
    {
        glm::vec3 dir = line.direction;
        glm::vec3 ori = line.origin;
        float A = glm::dot(dir, dir);
        float B = 2.f * (glm::dot(ori, dir) - glm::dot(dir, center));
        float C = glm::dot(center, center) + glm::dot(ori, ori) - 2.f * glm::dot(ori, center) - radius * radius;

        // discriminant
        float D = B * B - 4.f * A * C;
        if (D <= 0) return false;
        float t1 = (-B - sqrtf(D)) / (2.f * A);
        start = ori + t1 * dir;

     //   if (D == 0) return 1;

        float t2 = (-B + sqrtf(D)) / (2.f * A);
        end = ori + t2 * dir;

        if (getcolor) {
            glm::vec3 entry = start;
            if (t2 < t1) entry = end;
            glm::vec3 vector = glm::normalize(entry - center);
            color.x = 0.2f;
            color.y = acos(vector.z) / M_PI;
            color.z = atan(vector.y / vector.x) / (M_PI / 2.f);
        }

        return true;
    }

    void vaoGeneration(int stackCount, int sectorCount) {
        std::vector<Vertex> vertices;

        GLfloat x, y, z, xy;                              // vertex position
        GLfloat nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
       // GLfloat s, t;                                     // vertex texCoord

        GLfloat sectorStep = 2.f * glm::pi<float>() / (float)sectorCount;
        GLfloat stackStep = glm::pi<float>() / (float)stackCount;
        GLfloat sectorAngle, stackAngle;

        for (int i = 0; i <= stackCount; ++i)
        {
            stackAngle = glm::pi<float>() / 2 - i * stackStep;        // starting from pi/2 to -pi/2
            xy = radius * cosf(stackAngle);             // r * cos(u)
            z = radius * sinf(stackAngle);              // r * sin(u)

            // add (sectorCount+1) vertices per stack
            // the first and last vertices have same position and normal, but different tex coords
            for (int j = 0; j <= sectorCount; ++j)
            {
                sectorAngle = j * sectorStep;           // starting from 0 to 2pi

                // vertex position (x, y, z)
                x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
                y = xy * sinf(sectorAngle);            // r * cos(u) * sin(v)

                // normalized vertex normal (nx, ny, nz)
                nx = x * lengthInv;
                ny = y * lengthInv;
                nz = z * lengthInv;
                vertices.push_back(Vertex{ glm::vec3(x,y,z), glm::vec3(nx, ny, nz), glm::vec3(0), glm::vec3(0), glm::vec3(0) });

                //// vertex tex coord (s, t) range between [0, 1]
                //s = (float)j / sectorCount;
                //t = (float)i / stackCount;
                //texCoords.push_back(s);
                //texCoords.push_back(t);
            }
        }

        std::vector<GLuint> indices;
        GLuint k1, k2;
        for (int i = 0; i < stackCount; ++i)
        {
            k1 = i * (sectorCount + 1);     // beginning of current stack
            k2 = k1 + sectorCount + 1;      // beginning of next stack

            for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
            {
                // 2 triangles per sector excluding first and last stacks
                // k1 => k2 => k1+1
                if (i != 0)
                {
                    indices.push_back(k1);
                    indices.push_back(k2);
                    indices.push_back(k1 + 1);
                }

                // k1+1 => k2 => k2+1
                if (i != (stackCount - 1))
                {
                    indices.push_back(k1 + 1);
                    indices.push_back(k2);
                    indices.push_back(k2 + 1);
                }
            }
        }
        glGenVertexArrays(1, &vao); // make VAO
        glBindVertexArray(vao);

        // copy interleaved vertex data (V/N/T) to VBO
       // GLuint vbo[2];
        GLuint vbo;
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

        // copy index data to VBO
        GLuint ibo;
        glGenBuffers(1, &ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);

        // The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
        // Stride is the distance in bytes between vertices
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, pos)));
        glEnableVertexAttribArray(0);

        // The normals should be retrieved from the same Vertex Buffer Object (glBindBuffer is optional)
        // The offset is different and the data should go to input 1 instead of 0
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, normal)));
        glEnableVertexAttribArray(1);

        indSize = indices.size();

    }

};

#endif