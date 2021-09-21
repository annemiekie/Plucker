#ifndef CUBE_H
#define CUBE_H

#include <GL/glew.h>

#include "ray.h"
#include <glm/vec3.hpp>
#include <glm/glm.hpp>

#include "geoObject.h"

class Cube : public GeoObject {
private:
    glm::vec3 bounds[2] = { glm::vec3(0), glm::vec3(0) };
    std::vector<std::vector<Ray>> quadLines = std::vector<std::vector<Ray>>(6, std::vector<Ray>(4));
    std::vector<std::vector<glm::vec3>> cornerPoints = std::vector<std::vector<glm::vec3>>(6, std::vector<glm::vec3>(4));

public:

    glm::vec3 size;


    Cube() { };

    Cube(glm::vec3 center, float cubesize) : GeoObject(center) {
        size = glm::vec3(cubesize);
        setBounds(center - cubesize / 2.f, center + cubesize / 2.f);
    };

    Cube(glm::vec3 minvec, glm::vec3 maxvec) {
        setBounds(minvec, maxvec);
    };

    glm::vec3 getBounds(int minmax) {
        return bounds[minmax];
    }

    void setBounds(glm::vec3 minvec, glm::vec3 maxvec) {
        bounds[0] = minvec;
        bounds[1] = maxvec;
        center = (bounds[0] + bounds[1]) / 2.f;
        size = glm::abs(bounds[0] - bounds[1]);
        makeCubeSideLines();
    }

    int getSgn(glm::vec3 mainDir) {
        return glm::dot(mainDir, glm::vec3(1)) < 0 ? 1 : 0;
    }

    int getIndex(glm::vec3 mainDir) {
        int b = getSgn(mainDir);
        for (int i = 0; i < 3; i++) {
            if (mainDir[i] != 0) return 3 * b + i;
        }
        return -1;
    }

    // CHECK IF THIS WORKS
    std::vector<Ray> getCubeSideLines(glm::vec3 mainDir) {
        return quadLines[getIndex(mainDir)];
    }

    std::vector<glm::vec3> getCubeCornerPoints(glm::vec3 mainDir) {
        return cornerPoints[getIndex(mainDir)];

    }

    void makeCubeSideLines() {
        makeCubeSideLine(glm::vec3(1, 0, 0));
        makeCubeSideLine(glm::vec3(0, 1, 0));
        makeCubeSideLine(glm::vec3(0, 0, 1));
        makeCubeSideLine(glm::vec3(-1, 0, 0));
        makeCubeSideLine(glm::vec3(0, -1, 0));
        makeCubeSideLine(glm::vec3(0, 0, -1));
    }
    // bounds need to be the right way 'around'
    void makeCubeSideLine(glm::vec3 mainDir) {
        std::vector<glm::vec3> cPoints(4);
        int max = -1;
        int b = getSgn(mainDir);

        glm::vec3 bigbounds[2] = { center - 1.5f * size, center + 1.5f * size };// {bounds[0], bounds[1]};// 
    
        bool first = true;
        int diri = -1;
        for (int i = 0; i < 3; i++) {
            if (mainDir[i] == 0) {
                if (first) {
                    for (int j = 0; j < 4; j++) cPoints[j][i] = bigbounds[j % 3 == 0 ? 0 : 1][i];
                    first = false;
                }
                else for (int j = 0; j < 4; j++) cPoints[j][i] = bigbounds[j < 2 ? 0 : 1][i];
            }
            else for (int j = 0; j < 4; j++) {
                diri = i;
                cPoints[j][i] = bounds[b][i];
            }
        }
    
        std::vector<Ray> sidelines;
        for (int i = 0; i < 4; i++) sidelines.push_back( Ray(cPoints[i], cPoints[(i+1)%4]) );
        quadLines[b*3 + diri] = sidelines;
        cornerPoints[b * 3 + diri] = cPoints;
    }
    
    //// CHECK IF THIS WORKS
    bool intersectSide(glm::vec3 mainDir, const Ray& r, std::vector<Ray>& ignore, bool print) {
        int ind = getIndex(mainDir);
        int left = 0;
        int right = 0;
        for (int i = 0; i < 4; i++) {
            bool ign = false;
            for (auto& igray : ignore) {
                if (igray.equal(quadLines[ind][i], 1E-8)) {
                    ign = true; //-8
                    break;
                }
            }
            if (ign) continue;
            double sideVal = quadLines[ind][i].sideVal(r);
            if (sideVal < -1E-5) {
                if (print) std::cout << "Ray not in Box (wrong side of edge, should be >0): " << sideVal << std::endl <<std::endl;
                return false;
            }
            //if (abs(sideVal) < 1E-5) { left++; right++; }
            //else if (sideVal < 0) left++;
           // else right++;
        }
        return true;
       // if (right == 4) return true;// || left == 4
       // else {
       //     if (print)
       //     return false;
        //}
    
    }

    virtual bool intersect(const Ray& r, glm::vec3& start, glm::vec3& end, bool getcolor = false, glm::vec3& maindir = glm::vec3(0), glm::vec3& color = glm::vec3()) override  //const
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

        glm::vec3 tmin_xyz = glm::vec3(txmin, tymin, tzmin);
        glm::vec3 tmax_xyz = glm::vec3(txmax, tymax, tzmax);
        start = r.origin + (double)tmin * r.direction;
        end = r.origin + (double)tmax * r.direction;
        if (getcolor) {
            double x = glm::dot(maindir * maindir, tmin_xyz);
            glm::vec3 inter = r.origin + x * r.direction;
            glm::vec3 min = bounds[0] - size;
            float sz = size.x * 3.f;
            color = (inter - min) / sz;
            color.x = 0;
        }

        return true;
    };

    bool intersect(const Ray& r) const
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
        return true;
    };

    GLuint vaoGeneration() {
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

    GLuint vaoSideQuad(glm::vec3 mainDir) {
        int ind = getIndex(mainDir);
        std::vector<glm::vec3> lines;
        for (int i = 0; i < 4; i++) {
            lines.push_back(cornerPoints[ind][i]);
            lines.push_back(cornerPoints[ind][(i + 1) % 4]);
        }

        GLuint lineVAO, lineVBO;
        glGenVertexArrays(1, &lineVAO);
        glGenBuffers(1, &lineVBO);
        glBindVertexArray(lineVAO);
        glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*lines.size(), &lines[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);

        return lineVAO;
    }

};

#endif