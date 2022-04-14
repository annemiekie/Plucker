#ifndef CUBE_H
#define CUBE_H
#include <GL/glew.h>

#include "ray.h"
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include "axisAlignedSquare.h"
#include "vectorMath.h"
#include "geoObject.h"

class Cube : public GeoObject {
private:
    glm::dvec3 bounds[2] = { glm::dvec3(0), glm::dvec3(0) };
    glm::dvec3 bigbounds[2] = { glm::dvec3(0), glm::dvec3(0) };
    //std::vector<std::vector<Ray>> quadLines = std::vector<std::vector<Ray>>(6, std::vector<Ray>(4));
    //std::vector<std::vector<glm::vec3>> cornerPoints = std::vector<std::vector<glm::vec3>>(6, std::vector<glm::vec3>(4));
    std::vector<AxisAlignedSquare> sideSquares;

public:

    glm::dvec3 size;


    Cube() { };

    Cube(glm::dvec3 center, double cubesize) : GeoObject(center) {
        size = glm::dvec3(cubesize);
        setBounds(center - cubesize / 2.f, center + cubesize / 2.f);
    };

    Cube(glm::dvec3 minvec, glm::dvec3 maxvec) {
        setBounds(minvec, maxvec);
    };

    Cube(std::vector<glm::dvec3> points) {
        glm::dvec3 bbmin = glm::dvec3(1E10);
        glm::dvec3 bbmax = glm::dvec3(-1E10f);
        for (glm::dvec3& p : points) {
            for (int j = 0; j < 3; j++) {
                if (p[j] < bbmin[j]) bbmin[j] = p[j];
                if (p[j] > bbmax[j]) bbmax[j] = p[j];
            }
        }
        setBounds(bbmin, bbmax, false);
    };

    glm::dvec3 getBounds(int minmax) {
        return bounds[minmax];
    };

    //glm::dvec3 getBigBounds(int minmax) {
    //    return bigbounds[minmax];
    //}

    void setBounds(glm::vec3 minvec, glm::vec3 maxvec, bool makesquares = true) {
        bounds[0] = minvec;
        bounds[1] = maxvec;
        center = (bounds[0] + bounds[1]) / 2.;
        size = glm::abs(bounds[0] - bounds[1]);
        //vaoGeneration();
        if (makesquares) makeCubeSquares();
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

    AxisAlignedSquare getCubeSideSquare(glm::vec3 mainDir) {
        return sideSquares[getIndex(mainDir)];
    }

    //void makeCubeSquaresBig() {
    //    bigbounds[0] = center - 1.5 * size;
    //    bigbounds[1] = center + 1.5 * size;

    //    for (int j = 0; j < 2; j++) {
    //        for (int i = 0; i < 3; i++) {
    //            sideSquares.push_back(makeCubeSquareBig(i, j));
    //        }
    //    }
    //}

    void makeCubeSquares() {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 3; i++) {
                sideSquares.push_back(makeCubeSquare(i, j));
            }
        }
    }

    AxisAlignedSquare makeCubeSquare(int xyz, int pm) {
        
        glm::ivec3 normal;
        normal[xyz] = -1 + 2 * pm;

        glm::dvec3 point;
        point[xyz] = bounds[pm][xyz];

        glm::dvec2 minm, maxm;
        int c2 = 0;
        for (int c3 = 0; c3 < 3; c3++) {
            if (c3 != xyz) {
                minm[c2] = bounds[0][c3];
                maxm[c2] = bounds[1][c3];
                c2++;
            }
        }

        return AxisAlignedSquare(minm, maxm, point, normal);

        //glm::vec3 vmin, vside1, vside2, vopp;
        //vmin[xyz] = bounds[pm][xyz];
        //vmin[(xyz + 1) % 3] = bounds[0][(xyz + 1) % 3];
        //vmin[(xyz + 2) % 3] = bounds[0][(xyz + 2) % 3];

        //vside1[xyz] = vmin[xyz];
        //vside1[(xyz + 1) % 3] = bounds[1][(xyz + 1) % 3];
        //vside1[(xyz + 2) % 3] = bounds[0][(xyz + 2) % 3];

        //vside2[xyz] = vmin[xyz];
        //vside2[(xyz + 1) % 3] = bounds[0][(xyz + 1) % 3];
        //vside2[(xyz + 2) % 3] = bounds[1][(xyz + 2) % 3];

        //vopp = vmin;
        //vopp[xyz] = bounds[1 - pm][xyz];

        //return Square(vmin, vside1, vside2, vopp);
    }

    //// xyz = x, y or z direction, pm = plus or minus
    //Square makeCubeSquareBig(int xyz, int pm) {
    //    glm::vec3 vmin, vside1, vside2, vopp;
    //    vmin[xyz] = bounds[pm][xyz];
    //    vmin[(xyz + 1) % 3] = bigbounds[0][(xyz + 1) % 3];
    //    vmin[(xyz + 2) % 3] = bigbounds[0][(xyz + 2) % 3];

    //    vside1[xyz] = vmin[xyz];
    //    vside1[(xyz + 1) % 3] = bigbounds[1][(xyz + 1) % 3];
    //    vside1[(xyz + 2) % 3] = bigbounds[0][(xyz + 2) % 3];

    //    vside2[xyz] = vmin[xyz];
    //    vside2[(xyz + 1) % 3] = bigbounds[0][(xyz + 1) % 3];
    //    vside2[(xyz + 2) % 3] = bigbounds[1][(xyz + 2) % 3];

    //    vopp = vmin;
    //    vopp[xyz] = bounds[1 - pm][xyz];

    //    return Square(vmin, vside1, vside2, vopp);
    //}

    //// CHECK IF THIS WORKS
    bool intersectSide(glm::vec3 mainDir, const Ray& r) {
        if (glm::dot(glm::vec3(r.direction), mainDir) < 0) return false;
        return getCubeSideSquare(mainDir).inBounds(r, 1E-8);
    }

    glm::vec3 intersectSidePoint(glm::vec3 mainDir, const Ray& r) {
        return getCubeSideSquare(mainDir).rayIntersection(r);
        //int b = getSgn(mainDir);
        //double ttest = glm::dot(((bounds[b] - (glm::vec3)r.origin) / (glm::vec3)r.direction), glm::abs(mainDir));
        //glm::vec3 intersectpt = r.origin + r.direction * ttest;
        //return intersectpt;
    }

    bool intersectSegmSegm(glm::vec2& s1p1, glm::vec2& s1p2, glm::vec2& s2p1, glm::vec2& s2p2) {
        float d = (s1p2.x - s1p1.x) * (s2p2.y - s2p1.y) - (s1p2.y - s1p1.y) * (s2p2.x - s2p1.x);
        if (d == 0) return false;
        float q = (s1p1.y - s2p1.y) * (s2p2.x - s2p1.x) - (s1p1.x - s2p1.x) * (s2p2.y - s2p1.y);
        float r = q / d;
        q = (s1p1.y - s2p1.y) * (s1p2.x - s1p1.x) - (s1p1.x - s2p1.x) * (s1p2.y - s1p1.y);
        float s = q / d;
        if (r < 0 || r > 1 || s < 0 || s > 1) return false;
        return true;
    }


    //bool intersectSideSwathSimple(glm::vec3& vpos, glm::vec3& v1pos, glm::vec3& v2pos, glm::vec3& mainDir, std::vector<glm::vec3>& intersects = std::vector<glm::vec3>()) {
    //    glm::vec3 intr1 = intersectSidePoint(mainDir, Ray(vpos, v1pos));
    //    glm::vec3 intr2 = intersectSidePoint(mainDir, Ray(vpos, v2pos));
    //    intersects = { intr1, intr2 };
    //    glm::vec3 invMain = 1.f - glm::abs(mainDir);
    //    glm::vec3 sidemin = invMain * bigbounds[0];
    //    glm::vec3 sidemax = invMain * bigbounds[1];
    //    intr1 = intr1 * invMain;
    //    if (glm::all(glm::lessThanEqual(intr1, sidemax)) && glm::all(glm::greaterThanEqual(intr1, sidemin))) return true;
    //    intr2 = intr2 * invMain;
    //    if (glm::all(glm::lessThanEqual(intr2, sidemax)) && glm::all(glm::greaterThanEqual(intr2, sidemin))) return true;
    //    return false;
    //}

    // should do this with square
    bool intersectSideBB(std::vector<Ray>& rays, glm::vec3& mainDir) {
        std::vector<glm::dvec3> points;
        for (Ray& r : rays) points.push_back(intersectSidePoint(mainDir, r));
        Cube bb(points);
        for (int j = 0; j < 3; j++) {
            if (!mainDir[j]) {
                if (bb.bounds[1][j] < bounds[0][j] || bb.bounds[0][j] > bounds[1][j]) return false;
            }
        }
        return true;
    }

    bool intersectSideSwath(glm::dvec3& vpos, glm::dvec3& v1pos, glm::dvec3& v2pos, glm::vec3& mainDir, std::vector<glm::vec3>& intersects = std::vector<glm::vec3>()) {

        glm::dvec3 intr1 = intersectSidePoint(mainDir, Ray(vpos, v1pos));
        glm::dvec3 intr2 = intersectSidePoint(mainDir, Ray(vpos, v1pos));
        intersects = { intr1, intr2 };
        glm::dvec3 invMain = 1.f - glm::abs(mainDir);
        glm::dvec3 sidemin = invMain * bounds[0];
        glm::dvec3 sidemax = invMain * bounds[1];
        intr1 = intr1 * invMain;
        if (glm::all(glm::lessThanEqual(intr1, sidemax)) && glm::all(glm::greaterThanEqual(intr1, sidemin))) return true;
        intr1 = intr1 * invMain;
        if (glm::all(glm::lessThanEqual(intr2, sidemax)) && glm::all(glm::greaterThanEqual(intr2, sidemin))) return true;
        glm::vec2 int2r1, int2r2;
        glm::vec2 smin, smax;
        int count = 0;
        for (int i = 0; i < 3; i++) {
            if (invMain[i]) {
                if (intr1[i] < sidemin[i] && intr2[i] < sidemin[i]) return false;
                else if (intr1[i] > sidemax[i] && intr2[i] > sidemax[i]) return false;
                int2r1[count] = intr1[i];
                int2r2[count] = intr2[i];
                smin[count] = sidemin[i];
                smax[count] = sidemax[i];
                count++;
            }
        }
        if (intersectSegmSegm(int2r1, int2r2, smin, glm::vec2{ smin.x, smax.y })) return true;
        if (intersectSegmSegm(int2r1, int2r2, smin, glm::vec2{ smin.y, smax.x })) return true;
        if (intersectSegmSegm(int2r1, int2r2, smax, glm::vec2{ smin.x, smax.y })) return true;
        if (intersectSegmSegm(int2r1, int2r2, smax, glm::vec2{ smin.y, smax.x })) return true;

        return false;
    }

    virtual bool intersect(const Ray& r, double& tmin, double& tmax, bool getcolor = false, glm::vec3& maindir = glm::vec3(0), glm::vec3& color = glm::vec3()) override  //const
    {
        bool sign[3];
        sign[0] = (r.invdir.x < 0);
        sign[1] = (r.invdir.y < 0);
        sign[2] = (r.invdir.z < 0);
        double txmin, txmax, tymin, tymax, tzmin, tzmax;
        tmin = -INFINITY;
        tmax = INFINITY;

        if (r.direction.x == 0) {
            if (r.origin.x < bounds[sign[0]].x) return false;
            if (r.origin.x > bounds[1 - sign[0]].x) return false;
        }
        else {
            txmin = (bounds[sign[0]].x - r.origin.x) * r.invdir.x;
            txmax = (bounds[1 - sign[0]].x - r.origin.x) * r.invdir.x;
            tmin = txmin;
            tmax = txmax;
        }

        if (r.direction.y == 0) {
            if (r.origin.y < bounds[sign[0]].y) return false;
            if (r.origin.y > bounds[1 - sign[0]].y) return false;
        }
        else {
            tymin = (bounds[sign[1]].y - r.origin.y) * r.invdir.y;
            tymax = (bounds[1 - sign[1]].y - r.origin.y) * r.invdir.y;
            if ((tmin > tymax) || (tymin > tmax))
                return false;
            if (tymin > tmin)
                tmin = tymin;
            if (tymax < tmax)
                tmax = tymax;
        }
        if (r.direction.z == 0) {
            if (r.origin.z < bounds[sign[0]].z) return false;
            if (r.origin.z > bounds[1 - sign[0]].z) return false;
        }
        else {
            tzmin = (bounds[sign[2]].z - r.origin.z) * r.invdir.z;
            tzmax = (bounds[1 - sign[2]].z - r.origin.z) * r.invdir.z;

            if ((tmin > tzmax) || (tzmin > tmax))
                return false;
            if (tzmin > tmin)
                tmin = tzmin;
            if (tzmax < tmax)
                tmax = tzmax;
        }

        if (getcolor) {
            glm::dvec3 tmin_xyz = glm::dvec3(txmin, tymin, tzmin);
            glm::dvec3 tmax_xyz = glm::dvec3(txmax, tymax, tzmax);
            double x = dot_fd(maindir * maindir, tmin_xyz);
            glm::dvec3 inter = r.origin + x * r.direction;
            glm::dvec3 min = bounds[0] - size;
            double sz = size.x * 3.;
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
        float txmin, txmax, tymin, tymax, tzmin, tzmax;
        float tmin = -INFINITY;
        float tmax = INFINITY;

        if (r.direction.x == 0) {
            if (r.origin.x < bounds[sign[0]].x) return false;
            if (r.origin.x > bounds[1 - sign[0]].x) return false;
        }
        else {
            txmin = (bounds[sign[0]].x - r.origin.x) * r.invdir.x;
            txmax = (bounds[1 - sign[0]].x - r.origin.x) * r.invdir.x;
            tmin = txmin;
            tmax = txmax;
        }

        if (r.direction.y == 0) {
            if (r.origin.y < bounds[sign[0]].y) return false;
            if (r.origin.y > bounds[1 - sign[0]].y) return false;
        }
        else {
            tymin = (bounds[sign[1]].y - r.origin.y) * r.invdir.y;
            tymax = (bounds[1 - sign[1]].y - r.origin.y) * r.invdir.y;
            if ((tmin > tymax) || (tymin > tmax))
                return false;
            if (tymin > tmin)
                tmin = tymin;
            if (tymax < tmax)
                tmax = tymax;
        }
        if (r.direction.z == 0) {
            if (r.origin.z < bounds[sign[0]].z) return false;
            if (r.origin.z > bounds[1 - sign[0]].z) return false;
        }
        else {
            tzmin = (bounds[sign[2]].z - r.origin.z) * r.invdir.z;
            tzmax = (bounds[1 - sign[2]].z - r.origin.z) * r.invdir.z;

            if ((tmin > tzmax) || (tzmin > tmax))
                return false;
        }
        return true;
    };

    void vaoGeneration() {
        std::vector<GLdouble> lineSegments = {
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

        GLuint lineVBO;
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &lineVBO);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble) * lineSegments.size(), &lineSegments[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(GLdouble), (void*)0);
        //return lineVAO;
    };

    GLuint vaoSideQuad(glm::vec3 mainDir) {
        //int ind = getIndex(mainDir);
        Square sq = getCubeSideSquare(mainDir);
        std::vector<glm::vec3> lines;

        for (int i = 0; i < 4; i++) {
            lines.push_back(sq.vertices[i]);
            lines.push_back(sq.vertices[(i + 1) % 4]);
        }

        GLuint lineVAO, lineVBO;
        glGenVertexArrays(1, &lineVAO);
        glGenBuffers(1, &lineVBO);
        glBindVertexArray(lineVAO);
        glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * lines.size(), &lines[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);

        return lineVAO;
    };

    virtual void draw() override {
        glBindVertexArray(vao);
        glDrawArrays(GL_LINES, 0, 24);
    };

};

#endif