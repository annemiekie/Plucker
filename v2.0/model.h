#ifndef MODEL_H
#define MODEL_H

#include <GL/glew.h>

#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <embree3/rtcore.h>
#include <tiny_obj_loader.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "vertex.h"
#include "edge.h"
#include "cube.h"
#include "sphere.h"
#include <set>
#include <utility>
#include "primitive.h"
#include "RobinHood.h"
#include "tetrahedron.h"

class Model {
public:
	GLuint vao = 0;
	GLuint vbo = 0;
	GLuint vao2 = 0;

	//enlarged triangles
	std::vector<Vertex> verticesL = std::vector<Vertex>();

	std::vector<Vertex> vertices = std::vector<Vertex>();
	std::vector<glm::vec3> verticesIndexed = std::vector<glm::vec3>();

	std::vector<glm::uint> indices = std::vector<glm::uint>();
	std::set<Edge> edges = std::set<Edge>();

	std::set<Edge> silhouetteEdges = std::set<Edge>();

	std::vector<std::vector<Edge>> edgesPerTriangle = std::vector<std::vector<Edge>>();
	std::vector<std::vector<Edge>> edgesPerVertex = std::vector<std::vector<Edge>>();

	std::vector<std::vector<int>> triPerVertex;
	std::vector<glm::vec3> normalPerTri = std::vector<glm::vec3>();

	// need 6 of these
	bool cacheEE = false;
	bool cacheEEE = false;
	std::vector<robin_hood::unordered_map<uint64_t, bool>> edgeEdgeCombis;// = std::vector<std::unordered_map<long, bool>>();
	//std::vector<std::vector<uint64_t>> edgeEdgeCombis;
	int cachehitee = 0;
	int cachehiteee = 0;

	robin_hood::unordered_map<uint64_t, bool> edgeEdgeEdgeCombis;

	int primsize = 0;
	Cube boundingBox;
	Cube boundingCube;
	Sphere boundingSphere;
	glm::vec3 center = glm::vec3(0, 0, 0);
	float radius = 0.f;

	RTCDevice device;
	RTCScene scene;
	struct RTCIntersectContext context;

	Model() {};

	~Model() {
		rtcReleaseScene(scene);
	}

	Model(const char* filename, bool indexed = true, bool cacheEE = false, bool cacheEEE = false) :
		cacheEE(cacheEE), cacheEEE(cacheEEE) {
		loadModelFromFile(filename);
		//findPotentialSilhouettes(glm::vec3(1, 0, 0));
		setUpEmbreeTracer();
		addGeometry(indexed);
		createVAO();
		createIBO();
		//enlargeModel();
		if (cacheEE) makeCacheEE();

	}

	void makeCacheEE() {
		for (int i = 0; i < 6; i++) {
			edgeEdgeCombis.push_back(robin_hood::unordered_map<uint64_t, bool>());
		}
	}

	void enlargeModel() {
		//Model newmodel = Model(*this);
		verticesL = std::vector<Vertex>(vertices.size());
		for (int i = 0; i < vertices.size(); i++)
			verticesL[i].pos = vertices[i].pos + glm::normalize(vertices[i].pos - vertices[i].center) * 0.00001f;

		detachGeometry();
		rtcReleaseScene(scene);
		setUpEmbreeTracer();
		addGeometry(false);
	}

	void loadModelFromFile(const char* filename) {
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string err, warn;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename)) {
			std::cerr << err << std::endl;
			return;
		}

		int ind = 0;
		int count = 0;
		float n = 1000.f;
		float maxx = -n, maxy = -n, maxz = -n;
		float minx = n, miny = n, minz = n;

		verticesIndexed = std::vector<glm::vec3>(attrib.vertices.size() / 3);
		triPerVertex = std::vector<std::vector<int>>(attrib.vertices.size() / 3);
		edgesPerVertex = std::vector<std::vector<Edge>>(attrib.vertices.size() / 3);

		//for (int i = 0; i < attrib.vertices.size(); i+=3) {
		//	glm::vec3 vertex = {
		//		attrib.vertices[i],
		//		attrib.vertices[i + 1],
		//		attrib.vertices[i + 2]
		//	};
		//	vertices2[i/3] = vertex;
		//}

		std::vector<int> verticespertri(3);
		//glm::vec3 normalTri(0);
		std::vector<glm::vec3> centroid(3);
		for (const auto& shape : shapes) {
			primsize += shape.mesh.indices.size() / 3;

			for (const auto& index : shape.mesh.indices) {
				indices.push_back(index.vertex_index);
				Vertex vertex = {};

				// Retrieve coordinates for vertex by index
				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2]
				};

				vertex.bary = { count % 3 == 0, (count + 1) % 3 == 0, (count + 2) % 3 == 0 };

				if (vertex.pos.x > maxx) maxx = vertex.pos.x;
				else if (vertex.pos.x < minx) minx = vertex.pos.x;
				if (vertex.pos.y > maxy) maxy = vertex.pos.y;
				else if (vertex.pos.y < miny) miny = vertex.pos.y;
				if (vertex.pos.z > maxz) maxz = vertex.pos.z;
				else if (vertex.pos.z < minz) minz = vertex.pos.z;

				// Retrieve components of normal by index
				vertex.normal = {
					attrib.normals[3 * index.normal_index + 0],
					attrib.normals[3 * index.normal_index + 1],
					attrib.normals[3 * index.normal_index + 2]
				};

				vertex.color = { 1,1,1 };

				verticespertri[count % 3] = index.vertex_index;
				vertex.id = (1.f * ind);
				triPerVertex[index.vertex_index].push_back(ind - 1);
				verticesIndexed[index.vertex_index] = vertex.pos;
				vertices.push_back(vertex);

				if (count % 3 == 0) ind++;
				else if (count % 3 == 2) {
					normalPerTri.push_back(glm::normalize(glm::cross(vertices[vertices.size() - 3].pos - vertices[vertices.size() - 1].pos,
						vertices[vertices.size() - 2].pos - vertices[vertices.size() - 1].pos)));
					edgesPerTriangle.push_back(std::vector<Edge>(3));

					for (int x = 0; x < 3; x++) {
						int v1 = verticespertri[x];
						int v2 = verticespertri[(x + 1) % 3];
						if (v1 < v2) {
							int tmp = v1;
							v1 = v2;
							v2 = tmp;
						}
						Edge e = Edge({ { v1, v2 }, {}, -1 });
						std::pair<std::set<Edge>::iterator, bool> insert = edges.insert(e);
						(insert.first)->triangles.push_back(ind - 1);
						if (((insert.first)->index) == -1) {
							(insert.first)->index = edges.size() - 1;
							edgesPerVertex[verticespertri[x]].push_back(*(insert.first));
							edgesPerVertex[verticespertri[(x + 1) % 3]].push_back(*(insert.first));
						}
						edgesPerTriangle[edgesPerTriangle.size() - 1][x] = *(insert.first);
					}
				}
				count++;
			}
		}

		for (int i = 0; i < vertices.size(); i += 3) {
			glm::vec3 centroid = (vertices[i].pos + vertices[i + 1].pos + vertices[i + 2].pos) / 3.f;
			vertices[i].center = centroid;
			vertices[i + 1].center = centroid;
			vertices[i + 2].center = centroid;
		}

		boundingBox = Cube(glm::vec3(minx, miny, minz), glm::vec3(maxx, maxy, maxz));
		center = glm::vec3((minx + maxx) / 2.f, (miny + maxy) / 2.f, (minz + maxz) / 2.f);
		radius = glm::length(glm::vec3(maxx, maxy, maxz) - center);
		float maxside = std::max(std::max(maxx - minx, maxy - miny), maxz - minz);
		boundingCube = Cube(center, maxside);
		boundingSphere = Sphere(center, glm::length(boundingCube.getBounds(1) - center));
	};

	RTCGeometry makeEmbreeGeom(bool indexed = true) {
		RTCGeometry model_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
		int vertexsize = indexed ? verticesIndexed.size() : vertices.size();
		float* vertices_embr = (float*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), vertexsize);
		unsigned* triangles_embr = (unsigned*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), indices.size() / 3);
		for (int i = 0; i < vertexsize; i++) {
			vertices_embr[i * 3] = indexed ? verticesIndexed[i][0] : verticesL[i].pos.x;
			vertices_embr[i * 3 + 1] = indexed ? verticesIndexed[i][1] : verticesL[i].pos.y;
			vertices_embr[i * 3 + 2] = indexed ? verticesIndexed[i][2] : verticesL[i].pos.z;
		}
		for (int i = 0; i < indices.size(); i++)
			triangles_embr[i] = indexed ? indices[i] : i;
		return model_geometry;
	}

	void setUpEmbreeTracer() {
		device = rtcNewDevice(NULL);
		scene = rtcNewScene(device);
	}

	void detachGeometry() {
		rtcDetachGeometry(scene, 0);
	}

	void addGeometry(bool indexed = true) {
		RTCGeometry model_geometry = makeEmbreeGeom(indexed);
		rtcCommitGeometry(model_geometry);
		rtcAttachGeometry(scene, model_geometry);
		rtcReleaseGeometry(model_geometry);
		rtcCommitScene(scene);
		rtcInitIntersectContext(&context);
	}

	//void findPotentialSilhouettes(glm::vec3 maindir) {
	//	std::vector<glm::vec3> cpoints = boundingCube.getCubeCornerPoints(maindir);
	//	for (auto &e : edges) {
	//		bool check = false;
	//		if (e.triangles.size() == 1) check = true;
	//		else {
	//			int count = 0;
	//			for (auto &c : cpoints) {
	//				glm::vec3 halfway = 0.5f * (verticesIndexed[e.v[0]] + verticesIndexed[e.v[1]]);
	//				bool dot1 = glm::dot(glm::normalize(halfway - c), normalPerTri[e.triangles[0]]) < 0;
	//				bool dot2 = glm::dot(glm::normalize(halfway - c), normalPerTri[e.triangles[1]]) < 0;
	//				if (dot1 != dot2) {
	//					check = true;
	//					break;
	//				}
	//				else if (!dot1) count++;
	//			}
	//			if (!check && count > 0 && count < 4) check = true;
	//		}
	//		if (check) {
	//			silhouetteEdges.insert(e);
	//			for (auto t : e.triangles) {
	//				vertices[3 * t].selected = -1.f;
	//				vertices[3 * t + 1].selected = -1.f;
	//				vertices[3 * t + 2].selected = -1.f;
	//			}
	//		}
	//	}
	//}


	bool getIntersectionEmbree(const Ray& ray, int& primIndex, float& depth, bool backculling = false) {
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
			// check orientation
			if (backculling) {
				glm::vec3 norm = vertices[rayhit.hit.primID * 3].normal;
				if (glm::dot(norm, (glm::vec3)ray.direction) > 0) return false;
			}
			primIndex = rayhit.hit.primID;
			depth = rayhit.ray.tfar;
			return true;
		}
		//std::cout << "The ray did not collide with the object." << std::endl;

		return false;
	}

	std::vector<Ray> getEdgeRaysForPrim(int prim) {
		std::vector<Edge> triEdges = edgesPerTriangle[prim];
		std::vector<Ray> edgeRays;
		glm::vec3 normal = normalPerTri[prim];
		glm::vec3 center = vertices[3 * prim].center;
		Ray normalInvert = Ray(center + normal, center);
		for (int i = 0; i < 3; i++) {
			Ray triEdge = Ray(verticesIndexed[triEdges[i].v[0]], verticesIndexed[triEdges[i].v[1]], triEdges[i].index);
			if (triEdge.side(normalInvert)) triEdge.inverseDir();
			edgeRays.push_back(triEdge);
		}
		return edgeRays;
	}

	bool checkEdgeEdgeEdge(const Edge& e1, const Edge& e2, const Edge& e3) {
		uint64_t key;
		if (cacheEEE) {
			if (e1.index > e2.index) {
				if (e2.index > e3.index)		key = uint64_t(e1.index) << 42 | uint64_t(e2.index) << 21 | uint64_t(e3.index);
				else if (e1.index > e3.index)	key = uint64_t(e1.index) << 42 | uint64_t(e3.index) << 21 | uint64_t(e2.index);
				else							key = uint64_t(e3.index) << 42 | uint64_t(e1.index) << 21 | uint64_t(e2.index);
			}
			else {
				if (e3.index > e2.index)		key = uint64_t(e3.index) << 42 | uint64_t(e2.index) << 21 | uint64_t(e1.index);
				else if (e1.index > e3.index)	key = uint64_t(e2.index) << 42 | uint64_t(e1.index) << 21 | uint64_t(e3.index);
				else							key = uint64_t(e2.index) << 42 | uint64_t(e3.index) << 21 | uint64_t(e1.index);
			}
			if (edgeEdgeEdgeCombis.contains(key)) {
				cachehiteee++;
				return edgeEdgeEdgeCombis[key];
			}
		}

		std::vector<glm::vec3> e1v = { verticesIndexed[e1.v[0]], verticesIndexed[e1.v[1]] };
		std::vector<glm::vec3> e2v = { verticesIndexed[e2.v[0]], verticesIndexed[e2.v[1]] };
		Tetrahedron unboundTetra(e1v, e2v);

		//std::vector<float> d;
		//std::vector<glm::vec3> tri;
		//if (c0 != c[0] || c1 != c[1]) {
		//	std::vector<glm::vec3> e1v = { verticesIndexed[*e1.v[0]], verticesIndexed[*silhouetteEdges[c[0]].v[1]] };
		//	std::vector<glm::vec3> e2v = { model->verticesIndexed[*silhouetteEdges[c[1]].v[0]], model->verticesIndexed[*silhouetteEdges[c[1]].v[1]] };
		//	unboundTetra = Tetrahedron(e1v, e2v);
		//	c0 = c[0];
		//	c1 = c[1];
		//}

		std::vector<glm::vec3> e3v = { verticesIndexed[e3.v[0]], verticesIndexed[e3.v[1]] };
		if (unboundTetra.segmInTetra(e3v[0], e3v[1], 1E-7)) {
			if (cacheEEE) edgeEdgeEdgeCombis[key] = true;
			return true;
		}
		if (cacheEEE) edgeEdgeEdgeCombis[key] = false;
		return false;
	}

	bool checkEdgeEdgeCache(const Edge& e1, const Edge& e2, bool alldir, glm::vec3 maindir) {
		uint64_t key;
		int cacheIndex = 0;
		if (cacheEE) {
			if (!alldir) {
				int cacheIndex = maindir.x * 3 + maindir.y * 2 + maindir.z;
				cacheIndex < 0 ? cacheIndex += 6 : cacheIndex -= 3;
			}
			if (e1.index < e2.index) key = uint64_t(e1.index) << 32 | uint64_t(e2.index);
			else key = uint64_t(e2.index) << 32 | uint64_t(e1.index);
			if (edgeEdgeCombis[cacheIndex].contains(key)) {
				cachehitee++;
				return edgeEdgeCombis[cacheIndex].at(key);
			}
		}

		if (checkEdgeEdge(e1, e2, alldir, maindir)) {
			if (cacheEE) edgeEdgeCombis[cacheIndex][key] = true;
			return true;
		}
		if (cacheEE) edgeEdgeCombis[cacheIndex][key] = false;
		return false;
		//int8_t hit;
		//if (e1.index > e2.index) hit = edgeEdgeCombis[e1.index][e2.index];
		//else hit = edgeEdgeCombis[e2.index][e1.index];
		//if (hit >= 0) {
		//	cachehit++;
		//	return hit;
		//}
		//hit = checkEdgeEdge(e1, e2, alldir, maindir) ? 1 : 0;
		//if (e1.index > e2.index) edgeEdgeCombis[e1.index][e2.index] = hit;
		//else edgeEdgeCombis[e2.index][e1.index] = hit;
		//return hit == 1;
	}

	bool checkEdgeEdge(const Edge& e1, const Edge& e2, bool alldir, glm::vec3 maindir) {
		if (e1.v[0] == e2.v[0] || e1.v[1] == e2.v[1] ||
			e1.v[0] == e2.v[1] || e1.v[1] == e2.v[0])  return false;

		glm::vec3 e1v1 = verticesIndexed[e1.v[0]];
		glm::vec3 e1v2 = verticesIndexed[e1.v[1]];
		glm::vec3 e2v1 = verticesIndexed[e2.v[0]];
		glm::vec3 e2v2 = verticesIndexed[e2.v[1]];
		if (!alldir) {
			if (!boundingCube.intersectSideBB(std::vector<Ray> {Ray(e1v1, e2v1), Ray(e1v1, e2v2), Ray(e1v2, e2v1), Ray(e1v2, e2v2)}, maindir)) return false;
		}

		// edge 1 to 2
		bool side1;
		bool check12 = false;
		int silhEdgeCheck = checkSilhouetteEdge2(e1v1, e2, true, maindir, side1);
		if (silhEdgeCheck <= 0) {
			bool side2;
			silhEdgeCheck = checkSilhouetteEdge2(e1v2, e2, true, maindir, side2);
			if (silhEdgeCheck == 1);
			else if (silhEdgeCheck == -1) return false;
			else if (side1 == side2) return false;
		}

		// edge 2 to 1
		silhEdgeCheck = checkSilhouetteEdge2(e2v1, e1, true, maindir, side1);
		if (silhEdgeCheck <= 0) {
			bool side2;
			silhEdgeCheck = checkSilhouetteEdge2(e2v2, e1, true, maindir, side2);
			if (silhEdgeCheck == 1);
			else if (silhEdgeCheck == -1) return false;
			else if (side1 == side2) return false;
		}
		return true;
	}

	int checkSilhouetteEdge2(glm::vec3& vpos, const Edge& e, bool alldir, glm::vec3& maindir, bool& side) {

		if (e.triangles.size() == 2) {

			int v1 = e.v[0];
			int v2 = e.v[1];
			glm::vec3 v1pos = verticesIndexed[v1];
			glm::vec3 v2pos = verticesIndexed[v2];

			if (!alldir && (glm::dot(v1pos, maindir) > glm::dot(vpos, maindir) &&
				glm::dot(v2pos, maindir) > glm::dot(vpos, maindir))) return -1;

			// concave/convex check
			int t1 = e.triangles[0];
			int t2 = e.triangles[1];
			glm::vec3 vNotOnEdge1pos;
			glm::vec3 vNotOnEdge2pos;
			for (int i = 0; i < 3; i++) {
				int vNotOnEdge1 = indices[t1 * 3 + i];
				int vNotOnEdge2 = indices[t2 * 3 + i];
				if (vNotOnEdge1 != v1 && vNotOnEdge1 != v2) vNotOnEdge1pos = verticesIndexed[vNotOnEdge1];
				if (vNotOnEdge2 != v1 && vNotOnEdge2 != v2) vNotOnEdge2pos = verticesIndexed[vNotOnEdge2];
			}
			glm::vec3 normal1 = normalPerTri[e.triangles[0]];
			float concaveConvex = glm::dot(vNotOnEdge2pos - vNotOnEdge1pos, normal1);
			if (fabsf(concaveConvex) < 1E-8 || concaveConvex > 0) return -1; // in plane or concave

			// define plane of triangle 1 and 2
			glm::vec3 normal2 = normalPerTri[e.triangles[1]];
			float plane1Const = glm::dot(normal1, v1pos);
			float plane2Const = glm::dot(normal2, v2pos);
			bool pointToPlane1 = (glm::dot(normal1, vpos) - plane1Const) < 0;
			bool pointToPlane2 = (glm::dot(normal2, vpos) - plane2Const) < 0;
			if (pointToPlane1 != pointToPlane2) return 1;
			else { side = pointToPlane1; return 0; }

			// optional check: if the triangle projected via vertices of edge falls outside box bounds
		}
		else if (e.triangles.size() == 1) return 1;
		return -1;
	}

	bool checkSilhouetteEdge(glm::vec3& vpos, const Edge& e, bool alldir, glm::vec3& maindir) {


		if (e.triangles.size() == 2) {

			glm::vec3 v1 = verticesIndexed[e.v[0]];
			glm::vec3 v2 = verticesIndexed[e.v[1]];
			if (!alldir && (glm::dot(v1, maindir) > glm::dot(vpos, maindir) &&
				glm::dot(v2, maindir) > glm::dot(vpos, maindir))) return false;

			glm::vec3 halfway = 0.5f * (v1 + v2);
			float tri1normallinedot = glm::dot(glm::normalize(halfway - vpos), normalPerTri[e.triangles[0]]);
			float tri2normallinedot = glm::dot(glm::normalize(halfway - vpos), normalPerTri[e.triangles[1]]);

			// don't use edges in planes
			if (fabs(tri1normallinedot) < 1E-8 && fabs(tri2normallinedot) < 1E-8) return false;
			bool dot1 = tri1normallinedot < 0;
			bool dot2 = tri2normallinedot < 0;

			if (dot1 == dot2) return false;
			glm::vec3 jitter = glm::normalize(halfway - vpos - vertices[e.triangles[0] * 3].center);
			float dotjittertri1 = glm::dot(jitter, normalPerTri[e.triangles[0]]);
			float dotjittertri2 = glm::dot(jitter, normalPerTri[e.triangles[1]]);
			if (dotjittertri1 > 0 || dotjittertri2 > 0) {
				return true;
			}

		}
		else if (e.triangles.size() == 1) return true;
		return false;
	}

	//bool checkEdgeConnection(Edge& e1, Edge& e2) {
	//	if (e1.triangles.size() == 1) return true;
	//	else {

	//	}
	//}

	void findSilhouetteEdgesForTri(int prim, bool alldir, glm::vec3 maindir,
		std::vector<Edge>& silhouetteEdge,
		std::set<int>& checktris = std::set<int>()) {

		std::set<Edge> edgesToCheck = std::set<Edge>();

		if (checktris.size() != 0) 	for (int tri : checktris) for (Edge& edge : edgesPerTriangle[tri]) edgesToCheck.insert(edge);
		else edgesToCheck = edges;

		glm::vec3 invMaindir = glm::vec3(1.f) - maindir;
		glm::vec3 maindirMin = boundingBox.getBounds(0);
		glm::vec3 maindirMax = boundingBox.getBounds(1);

		for (const Edge& e : edgesToCheck) {
			// If edge is not in potential silhouettes
			//if (!alldir && silhouetteEdges.find(e) == silhouetteEdges.end()) continue;

			// Don't use same triangle as silhouette edge.
			if (e.triangles[0] == prim) continue;
			if (e.triangles.size() == 2) if (e.triangles[1] == prim) continue;

			// check if edge lies on correct side of triangle
			if (glm::dot(normalPerTri[prim], verticesIndexed[e.v[0]] - vertices[3 * prim].center) < 0 &&
				glm::dot(normalPerTri[prim], verticesIndexed[e.v[1]] - vertices[3 * prim].center) < 0) continue;

			glm::vec3 v1pos = verticesIndexed[e.v[0]];
			glm::vec3 v2pos = verticesIndexed[e.v[1]];

			bool sideCheck = false;
			bool intersect;
			std::vector<glm::vec3> intersectpts;
			bool toAdd = false;
			for (int i = 0; i < 3; i++) {
				bool side = false;
				int checkEdge = checkSilhouetteEdge2(vertices[3 * prim + i].pos, e, alldir, maindir, side);
				if (!alldir) {
					Ray r1 = Ray(v1pos, vertices[3 * prim + i].pos);
					Ray r2 = Ray(v2pos, vertices[3 * prim + i].pos);

					intersect = boundingCube.intersectSide(maindir, r1) || boundingCube.intersectSide(maindir, r2);
					intersectpts.push_back(boundingCube.intersectSidePoint(maindir, r1));
					intersectpts.push_back(boundingCube.intersectSidePoint(maindir, r2));
				}
				if (checkEdge == 0 && (i == 0 || side == sideCheck)) {
					sideCheck = side;
					continue;
				}
				else if (checkEdge >= 0) {
					if (alldir) {
						silhouetteEdge.push_back(e);
						break;
					}
					else {
						toAdd = !intersect;
						if (intersect) {
							silhouetteEdge.push_back(e);
							break;
						}
					}
				}
			}
			if (toAdd) {
				Cube bb(intersectpts);
				for (int i = 0; i < 3; i++) {
					if (invMaindir[i] == 1) {
						if (bb.getBounds(1)[i] < maindirMin[i] || bb.getBounds(0)[i] > maindirMax[i]) {
							toAdd = false;
							break;
						}
					}
				}
				if (toAdd) silhouetteEdge.push_back(e);
				//glm::vec2 bbmin2;
				//glm::vec2 bbmax2;
				//int count = 0;
				//for (int i = 0; i < 3; i++) {
				//	if (invMaindir[i] == 1) {
				//		bbmin2[count] = bbmin[i];
				//		bbmax2[count] = bbmax[i];
				//		count++;
				//	}
				//}
				//if (bbmin2.x > maindirMin.x && bbmin2.x < maindirMax.x) 
			}

			// if no swath intersect check bounding box
		}
	}

	double getIntersectionDepthForPrim(const int i, const Ray& r) {
		glm::dvec3 v0 = verticesL[i].pos;
		glm::dvec3 v1 = verticesL[i + 1].pos;
		glm::dvec3 v2 = verticesL[i + 2].pos;

		glm::dvec3 N = glm::normalize(glm::cross(v2 - v0, v1 - v0));
		//glm::dvec3 L = glm::vec3(-3, 3, -3);
		//diff = std::max(glm::dot(-glm::normalize(N), glm::normalize(L)), 0.);
		double ndotdir = glm::dot(N, r.direction);
		double t = glm::dot(N, v0 - r.origin) / ndotdir;
		return t;
	}

	bool getIntersectionWithPrim(int i, const Ray& r, float& depth) {
		Ray edge;
		// clockwise triangle edges
		edge = Ray(vertices[i].pos, vertices[i + 1].pos);
		if (edge.side(r)) return false;
		edge = Ray(vertices[i + 1].pos, vertices[i + 2].pos);
		if (edge.side(r)) return false;
		edge = Ray(vertices[i + 2].pos, vertices[i].pos);
		if (edge.side(r)) return false;
		depth = getIntersectionDepthForPrim(i, r);
		return true;
	};

	bool getIntersectionNoAcceleration(const Ray& r, int& primIndex, float& depth) {
		primIndex = -1;
		depth = 10000.f;
		for (int index = 0; index < primsize; index++) {
			float newdepth;
			if (getIntersectionWithPrim(index * 3, r, newdepth)) {
				if (newdepth < depth) {
					depth = newdepth;
					primIndex = index;
				}
			}
		}
		return primIndex >= 0;
	};


	void createIBO() {
		glGenVertexArrays(1, &vao2); // make VAO
		glBindVertexArray(vao2);

		// copy interleaved vertex data (V/N/T) to VBO
	   // GLuint vbo[2];
		GLuint vbo;
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, verticesIndexed.size() * sizeof(glm::vec3), &verticesIndexed[0], GL_STATIC_DRAW);

		// copy index data to VBO
		GLuint ibo;
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);

		// The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
		// Stride is the distance in bytes between vertices
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
		glEnableVertexAttribArray(0);

	}

	void draw() {
		glBindVertexArray(vao);
		glDrawArrays(GL_TRIANGLES, 0, vertices.size());
	}

	void clearColors(glm::vec3& color) {
		//for (Vertex& v : vertices) v.color = color;
		for (int i = 0; i < vertices.size(); i++) vertices[i].color = color;
	}

	template <class T>
	void changePrimColors(T& setOfprims, glm::vec3& color) {
		for (auto i : setOfprims) changePrimColor(i, color);
	}

	void changePrimColor(int prim, glm::vec3& color) {
		for (int j = 0; j < 3; j++) {
			vertices[3 * prim + j].color = color;
		}
	}
	void changeSelected() {
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);
		//glBindVertexArray(vao);
		//glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, selected)));
		//glEnableVertexAttribArray(4);
	}

	void createVAO() {
		//////////////////// Create Vertex Buffer Object
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

		// Bind vertex data to shader inputs using their index (location)
		// These bindings are stored in the Vertex Array Object
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		// The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
		// Stride is the distance in bytes between vertices
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, pos)));
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, normal)));
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, bary)));
		glEnableVertexAttribArray(2);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, center)));
		glEnableVertexAttribArray(3);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, color)));
		glEnableVertexAttribArray(4);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, id)));
		glEnableVertexAttribArray(5);
	}

};

#endif