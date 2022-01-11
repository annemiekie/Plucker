#ifndef CUBE_H
#define CUBE_H
#include <GL/glew.h>

#include "ray.h"
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include "square.h"

#include "geoObject.h"

class Cube : public GeoObject {
private:
    glm::vec3 bounds[2] = { glm::vec3(0), glm::vec3(0) };
   // glm::vec3 bigbounds[2] = { glm::vec3(0), glm::vec3(0) };
    //std::vector<std::vector<Ray>> quadLines = std::vector<std::vector<Ray>>(6, std::vector<Ray>(4));
    //std::vector<std::vector<glm::vec3>> cornerPoints = std::vector<std::vector<glm::vec3>>(6, std::vector<glm::vec3>(4));
    std::vector<Square> sideSquares;

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

    Cube(std::vector<glm::vec3> points) {
        glm::vec3 bbmin = glm::vec3(1E10);
        glm::vec3 bbmax = glm::vec3(-1E10f);
        for (glm::vec3& p : points) {
            for (int j = 0; j < 3; j++) {
                if (p[j] < bbmin[j]) bbmin[j] = p[j];
                if (p[j] > bbmax[j]) bbmax[j] = p[j];
            }
        }
        setBounds(bbmin, bbmax, false);
    };

    glm::vec3 getBounds(int minmax) {
        return bounds[minmax];
    };

    //glm::vec3 getBigBounds(int minmax) {
    //    return bigbounds[minmax];
    //}

    void setBounds(glm::vec3 minvec, glm::vec3 maxvec, bool makesquares = true) {
        bounds[0] = minvec;
        bounds[1] = maxvec;
        center = (bounds[0] + bounds[1]) / 2.f;
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

    Square getCubeSideSquare(glm::vec3 mainDir) {
        return sideSquares[getIndex(mainDir)];
    }
    //std::vector<Ray> getCubeSideLines(glm::vec3 mainDir) {
    //    return quadLines[getIndex(mainDir)];
    //}

    //std::vector<glm::vec3> getCubeCornerPoints(glm::vec3 mainDir) {
    //    return cornerPoints[getIndex(mainDir)];

    //}

    void makeCubeSquares() {
        //bigbounds[0] = center - 1.5f * size;
        //bigbounds[1] = center + 1.5f * size;

        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 3; i++) {
                sideSquares.push_back(makeCubeSquare(i, j));
            }
        }
    }

    Square makeCubeSquare(int xyz, int pm) {
        glm::vec3 vmin, vside1, vside2, vopp;
        vmin[xyz] = bounds[pm][xyz];
        vmin[(xyz + 1) % 3] = bounds[0][(xyz + 1) % 3];
        vmin[(xyz + 2) % 3] = bounds[0][(xyz + 2) % 3];

        vside1[xyz] = vmin[xyz];
        vside1[(xyz + 1) % 3] = bounds[1][(xyz + 1) % 3];
        vside1[(xyz + 2) % 3] = bounds[0][(xyz + 2) % 3];

        vside2[xyz] = vmin[xyz];
        vside2[(xyz + 1) % 3] = bounds[0][(xyz + 1) % 3];
        vside2[(xyz + 2) % 3] = bounds[1][(xyz + 2) % 3];

        vopp = vmin;
        vopp[xyz] = bounds[1 - pm][xyz];

        return Square(vmin, vside1, vside2, vopp);
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


    // bounds need to be the right way 'around'
/*    void makeCubeSidelines(glm::vec3 mainDir) {
        std::vector<glm::vec3> cPoints(4);
        int max = -1;
        int b = getSgn(mainDir);

 
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
    }*/
    
    //// CHECK IF THIS WORKS
    bool intersectSide(glm::vec3 mainDir, const Ray& r) {
        if (glm::dot(glm::vec3(r.direction), mainDir) < 0) return false;
        return getCubeSideSquare(mainDir).inBounds(r, 1E-8);
        //int ind = getIndex(mainDir);
        //for (int i = 0; i < 4; i++) {
        //    bool ign = false;
        //     for (int j = 0; j < rayIgnoresize; j++) {
        //   //     if (lines[j].equal(quadLines[ind][i], 1E-8)) {
        //   // for (auto& igray : lines) {
        //   //      if (igray.equal(quadLines[ind][i], 1E-8)) {
        //            ign = true; //-8
        //            break;
        //        }
        //    }
        //    if (ign) continue;
        //    double sideVal = quadLines[ind][i].sideVal(r);
        //    if (sideVal < -1E-5) {
        //        if (print) std::cout << "Ray not in Box (wrong side of edge, should be >0): " << sideVal << std::endl <<std::endl;
        //        return false;
        //    }
        //}
        //return true;    
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
        q = (s1p1.y - s2p1.y) * (s1p2.x - s1p1.x) - (s1p1.x- s2p1.x) * (s1p2.y - s1p1.y);
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
        std::vector<glm::vec3> points;
        for (Ray& r : rays) points.push_back(intersectSidePoint(mainDir, r));
        Cube bb(points);
        for (int j = 0; j < 3; j++) {
            if (!mainDir[j]) {
                if (bb.bounds[1][j] < bounds[0][j] || bb.bounds[0][j] > bounds[1][j]) return false;
            }
        }
        return true;
    }

    bool intersectSideSwath(glm::vec3& vpos, glm::vec3& v1pos, glm::vec3& v2pos, glm::vec3& mainDir, std::vector<glm::vec3>& intersects = std::vector<glm::vec3>()) {

        glm::vec3 intr1 = intersectSidePoint(mainDir,Ray(vpos, v1pos));
        glm::vec3 intr2 = intersectSidePoint(mainDir,Ray(vpos, v1pos));
        intersects = { intr1, intr2 };
        glm::vec3 invMain = 1.f - glm::abs(mainDir);
        glm::vec3 sidemin = invMain * bounds[0];
        glm::vec3 sidemax = invMain * bounds[1];
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
                else if (intr1[i] > sidemax [i] && intr2[i] > sidemax[i]) return false;
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

    void vaoGeneration() {
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

        GLuint lineVBO;
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &lineVBO);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * lineSegments.size(), &lineSegments[0], GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
        //return lineVAO;
    };

    GLuint vaoSideQuad(glm::vec3 mainDir) {
        //int ind = getIndex(mainDir);
        Square sq = getCubeSideSquare(mainDir);
        std::vector<glm::vec3> lines;

        for (int i = 0; i < 4; i++) {
            lines.push_back(sq.cornerPoints[i]);
            lines.push_back(sq.cornerPoints[(i + 1) % 4]);
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

};

#endif