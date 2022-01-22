//#ifndef EMBREE_H
//#define EMBREE_H
//#include "raytracer.h"
//#include "model.h"
//
//class EmbreeTracer : public RayTracer {
//private:
//	RTCGeometry makeEmbreeGeom(Model* model, bool indexed = true) {
//		RTCGeometry model_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
//		int vertexsize = indexed ? model->verticesIndexed.size() : model->vertices.size();
//		float* vertices_embr = (float*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), vertexsize);
//		unsigned* triangles_embr = (unsigned*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), model->indices.size() / 3);
//		for (int i = 0; i < vertexsize; i++) {
//			vertices_embr[i * 3] = indexed ? model->verticesIndexed[i][0] : model->verticesL[i].pos.x;
//			vertices_embr[i * 3 + 1] = indexed ? model->verticesIndexed[i][1] : model->verticesL[i].pos.y;
//			vertices_embr[i * 3 + 2] = indexed ? model->verticesIndexed[i][2] : model->verticesL[i].pos.z;
//		}
//		for (int i = 0; i < model->indices.size(); i++)
//			triangles_embr[i] = indexed ? model->indices[i] : i;
//		return model_geometry;
//	}
//
//	void setUpEmbreeTracer() {
//		device = rtcNewDevice(NULL);
//		scene = rtcNewScene(device);
//	}
//
//public:
//	RTCDevice device;
//	RTCScene scene;
//	struct RTCIntersectContext context;
//
//	EmbreeTracer(Model* model, bool indexed = true) {
//		attachGeometry(model, indexed);
//	}
//
//	void detachGeometry() {
//		rtcDetachGeometry(scene, 0);
//	}
//
//	void attachGeometry(Model* model, bool indexed = true) {
//		RTCGeometry model_geometry = makeEmbreeGeom(model, indexed);
//		rtcCommitGeometry(model_geometry);
//		rtcAttachGeometry(scene, model_geometry);
//		rtcReleaseGeometry(model_geometry);
//		rtcCommitScene(scene);
//		rtcInitIntersectContext(&context);
//	}
//
//	//bool intersect()
//};
//#endif