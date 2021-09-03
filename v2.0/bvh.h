#pragma once
#include "accelerationStructure.h"
#include <math.h>
#include "ray.h"
#include "cube.h"
#include "vertex.h"

class BVH : public AccelerationStructure
{
    
public:
    BVH(const Cube boundingCube, const std::vector<Vertex>& vertices) {

    };
    //const Object* intersect(const Ray& ray, IsectData& isectData) const;
    ~BVH();
};