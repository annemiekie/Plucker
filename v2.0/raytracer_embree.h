#pragma once
#include "model.h"
#include <embree3/rtcore.h>
#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "embree3/common/math/bbox.h"

bool buildProgress(void* userPtr, double f) {
    return true;
}
void splitPrimitive(const RTCBuildPrimitive* prim, unsigned int dim, float pos, RTCBounds* lprim, RTCBounds* rprim, void* userPtr)
{
    assert(dim < 3);
    assert(prim->geomID == 0);
    *(embree::BBox3fa*)lprim = *(embree::BBox3fa*)prim;
    *(embree::BBox3fa*)rprim = *(embree::BBox3fa*)prim;
    (&lprim->upper_x)[dim] = pos;
    (&rprim->lower_x)[dim] = pos;
};

    class Raytracer_Embree {
    public:
        RTCDevice device;
        RTCScene scene;
        struct RTCIntersectContext context;


        Raytracer_Embree(Model* model) {
            device = rtcNewDevice(NULL);
            scene = rtcNewScene(device);
            
            RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
            float* vertices = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), model->vertices2.size());
            unsigned* triangles = (unsigned*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), model->indices.size() / 3);
            for (int i = 0; i < model->vertices2.size(); i++) {
                vertices[i * 3] = model->vertices2[i][0];
                vertices[i * 3 + 1] = model->vertices2[i][1];
                vertices[i * 3 + 2] = model->vertices2[i][2];
            }
            for (int i = 0; i < model->indices.size(); i++)
                triangles[i] = model->indices[i];

            rtcCommitGeometry(geom);
            rtcAttachGeometry(scene, geom);
            rtcReleaseGeometry(geom);
            rtcCommitScene(scene);
            rtcInitIntersectContext(&context);

           // std::vector<RTCBuildPrimitive> prims;
            //prim.
           // build(RTC_BUILD_QUALITY_MEDIUM, prims);
        }

        ~Raytracer_Embree() {
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        };

        /* prints the bvh4.triangle4v data structure */
        //void print_bvh4_triangle4v(BVH4::NodeRef node, size_t depth)
        //{
        //    if (node.isAABBNode())
        //    {
        //        BVH4::AABBNode* n = node.getAABBNode();
        //
        //        std::cout << "AABBNode {" << std::endl;
        //        for (size_t i = 0; i < 4; i++)
        //        {
        //            for (size_t k = 0; k < depth; k++) std::cout << "  ";
        //            std::cout << "  bounds" << i << " = " << n->bounds(i) << std::endl;
        //        }
        //
        //        for (size_t i = 0; i < 4; i++)
        //        {
        //            if (n->child(i) == BVH4::emptyNode)
        //                continue;
        //
        //            for (size_t k = 0; k < depth; k++) std::cout << "  ";
        //            std::cout << "  child" << i << " = ";
        //            print_bvh4_triangle4v(n->child(i), depth + 1);
        //        }
        //        for (size_t k = 0; k < depth; k++) std::cout << "  ";
        //        std::cout << "}" << std::endl;
        //    }
        //    else
        //    {
        //        size_t num;
        //        const Triangle4v* tri = (const Triangle4v*)node.leaf(num);
        //
        //        std::cout << "Leaf {" << std::endl;
        //        for (size_t i = 0; i < num; i++) {
        //            for (size_t j = 0; j < tri[i].size(); j++) {
        //                for (size_t k = 0; k < depth; k++) std::cout << "  ";
        //                std::cout << "  Triangle { v0 = (" << tri[i].v0.x[j] << ", " << tri[i].v0.y[j] << ", " << tri[i].v0.z[j] << "),  "
        //                    "v1 = (" << tri[i].v1.x[j] << ", " << tri[i].v1.y[j] << ", " << tri[i].v1.z[j] << "), "
        //                    "v2 = (" << tri[i].v2.x[j] << ", " << tri[i].v2.y[j] << ", " << tri[i].v2.z[j] << "), "
        //                    "geomID = " << tri[i].geomID(j) << ", primID = " << tri[i].primID(j) << " }" << std::endl;
        //            }
        //        }
        //        for (size_t k = 0; k < depth; k++) std::cout << "  ";
        //        std::cout << "}" << std::endl;
        //    }
        //}
//
        //void print_bvh(RTCScene scene)
        //{
        //    BVH4* bvh4 = nullptr;
        //
        //    /* if the scene contains only triangles, the BVH4 acceleration structure can be obtained this way */
        //    AccelData* accel = ((Accel*)scene)->intersectors.ptr;
        //    if (accel->type == AccelData::TY_BVH4)
        //        bvh4 = (BVH4*)accel;
        //    if (bvh4 == nullptr)
        //        throw std::runtime_error("cannot access BVH4 acceleration structure"); // will not happen if you use this Embree version
        //
        //      /* now lets print the entire hierarchy */
        //    print_bvh4_triangle4v(bvh4->root, 0);
        //}


        int intersect(Ray& ray) {
            struct RTCRayHit rayhit;
            rayhit.ray.org_x = ray.origin.x;
            rayhit.ray.org_y = ray.origin.y;
            rayhit.ray.org_z = ray.origin.z;
            rayhit.ray.dir_x = ray.direction.x;
            rayhit.ray.dir_y = ray.direction.y;
            rayhit.ray.dir_z = ray.direction.z;
            rayhit.ray.tnear = 0;
            rayhit.ray.tfar = std::numeric_limits<float>::infinity();
            rayhit.ray.mask = 0;
            rayhit.ray.flags = 0;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;
            rtcIntersect1(scene, &context, &rayhit);
            if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
            {
                //std::cout << "The ray collided with the geometry. The primitive ID is: " << rayhit.hit.primID << "Launch collision distance=" << rayhit.ray.tfar << std::endl;
                return rayhit.hit.primID;
            }
                //std::cout << "The ray did not collide with the object." << std::endl;

            return -1;
        }



   // private:

        //typedef void (*RTCSplitPrimitiveFunction) (const struct RTCBuildPrimitive* primitive, unsigned int dimension, float position, struct RTCBounds* leftBounds, struct RTCBounds* rightBounds, void* userPtr);



        struct NodeX
        {
            virtual float sah() = 0;
        };

        struct InnerNode : public NodeX
        {
            embree::BBox3fa bounds[2];
            NodeX* children[2];

            InnerNode() {
                bounds[0] = bounds[1] = embree::empty;
                children[0] = children[1] = nullptr;
            }

            float sah() {
                return 1.0f + (embree::area(bounds[0]) * children[0]->sah() + embree::area(bounds[1]) * children[1]->sah()) / embree::area(embree::merge(bounds[0], bounds[1]));
            }

            static void* create(RTCThreadLocalAllocator alloc, unsigned int numChildren, void* userPtr)
            {
                assert(numChildren == 2);
                void* ptr = rtcThreadLocalAlloc(alloc, sizeof(InnerNode), 16);
                return (void*) new (ptr) InnerNode;
            }

            static void  setChildren(void* nodePtr, void** childPtr, unsigned int numChildren, void* userPtr)
            {
                assert(numChildren == 2);
                for (size_t i = 0; i < 2; i++)
                    ((InnerNode*)nodePtr)->children[i] = (NodeX*)childPtr[i];
            }

            static void  setBounds(void* nodePtr, const RTCBounds** bounds, unsigned int numChildren, void* userPtr)
            {
                assert(numChildren == 2);
                for (size_t i = 0; i < 2; i++)
                    ((InnerNode*)nodePtr)->bounds[i] = *(const embree::BBox3fa*)bounds[i];
            }
        };

        struct LeafNode : public NodeX
        {
            unsigned id;
            embree::BBox3fa bounds;

            LeafNode(unsigned id, const embree::BBox3fa& bounds)
                : id(id), bounds(bounds) {}

            float sah() {
                return 1.0f;
            }

            static void* create(RTCThreadLocalAllocator alloc, const RTCBuildPrimitive* prims, size_t numPrims, void* userPtr)
            {
                assert(numPrims == 1);
                void* ptr = rtcThreadLocalAlloc(alloc, sizeof(LeafNode), 16);
                return (void*) new (ptr) LeafNode(prims->primID, *(embree::BBox3fa*)prims);
            }
        };

        void build(RTCBuildQuality quality, std::vector<RTCBuildPrimitive>& prims_i, size_t extraSpace = 0)
        {
            RTCBVH bvh = rtcNewBVH(device);

            std::vector<RTCBuildPrimitive> prims;
            prims.reserve(prims_i.size() + extraSpace);
            prims.resize(prims_i.size());

            /* settings for BVH build */
            RTCBuildArguments arguments = rtcDefaultBuildArguments();
            arguments.byteSize = sizeof(arguments);
            arguments.buildFlags = RTC_BUILD_FLAG_DYNAMIC;
            arguments.buildQuality = quality;
            arguments.maxBranchingFactor = 2;
            arguments.maxDepth = 1024;
            arguments.sahBlockSize = 1;
            arguments.minLeafSize = 1;
            arguments.maxLeafSize = 1;
            arguments.traversalCost = 1.0f;
            arguments.intersectionCost = 1.0f;
            arguments.bvh = bvh;
            arguments.primitives = prims.data();
            arguments.primitiveCount = prims.size();
            arguments.primitiveArrayCapacity = prims.capacity();
            arguments.createNode = InnerNode::create;
            arguments.setNodeChildren = InnerNode::setChildren;
            arguments.setNodeBounds = InnerNode::setBounds;
            arguments.createLeaf = LeafNode::create;
            arguments.splitPrimitive = splitPrimitive;
            arguments.buildProgress = buildProgress;
            arguments.userPtr = nullptr;

            /* we recreate the prims array here, as the builders modify this array */
            for (size_t j = 0; j < prims.size(); j++) prims[j] = prims_i[j];
            //           double t0 = getSeconds();
            NodeX* root = (NodeX*)rtcBuildBVH(&arguments);
           
            //           double t1 = getSeconds();
            //           const float sah = root ? root->sah() : 0.0f;
            //           std::cout << 1000.0f * (t1 - t0) << "ms, " << 1E-6 * double(prims.size()) / (t1 - t0) << " Mprims/s, sah = " << sah << " [DONE]" << std::endl;
            //       }

            //       rtcReleaseBVH(bvh);
            //   }
    
        }
    };
