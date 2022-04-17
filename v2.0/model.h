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
	std::vector<VertexVis> verticesVis;
	std::vector<glm::uint> indices;

	std::vector<Vertex*> vertices;
	std::vector<Primitive*> triangles;
	std::vector<Edge*> edgeVector;
	robin_hood::unordered_map <uint64_t, Edge*> edges;

	// need 6 of these
	bool cacheEE = false;
	bool cacheEEE = false;
	std::vector<robin_hood::unordered_map<uint64_t, bool>> edgeEdgeCombis;
	robin_hood::unordered_map<uint64_t, bool> edgeEdgeEdgeCombis;
	int cachehitee = 0;
	int cachehiteee = 0;

	int primsize = 0;
	Cube boundingBox;
	Cube boundingCube;
	Sphere boundingSphere;
	glm::dvec3 center = glm::dvec3(0);
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
		constructModelBounds();
		setUpEmbreeTracer();
		addGeometry(indexed);
		createVAO();
		//createIBO();
		//enlargeModel();
		if (cacheEE) makeCacheEE();

	}

	void makeCacheEE() {
		for (int i = 0; i < 6; i++) {
			edgeEdgeCombis.push_back(robin_hood::unordered_map<uint64_t, bool>());
		}
	}

	void constructModelBounds() {
		double maxx = -INFINITY, maxy = -INFINITY, maxz = -INFINITY;
		double minx = INFINITY, miny = INFINITY, minz = INFINITY;

		for (Vertex* vertex : vertices) {
			if (vertex->pos.x > maxx) maxx = vertex->pos.x;
			else if (vertex->pos.x < minx) minx = vertex->pos.x;
			if (vertex->pos.y > maxy) maxy = vertex->pos.y;
			else if (vertex->pos.y < miny) miny = vertex->pos.y;
			if (vertex->pos.z > maxz) maxz = vertex->pos.z;
			else if (vertex->pos.z < minz) minz = vertex->pos.z;
		}

		boundingBox = Cube(glm::dvec3(minx, miny, minz), glm::dvec3(maxx, maxy, maxz));
		center = glm::dvec3((minx + maxx) / 2., (miny + maxy) / 2., (minz + maxz) / 2.);
		radius = glm::length(glm::dvec3(maxx, maxy, maxz) - center);
		float maxside = std::max(std::max(maxx - minx, maxy - miny), maxz - minz);
		boundingCube = Cube(center, maxside);
		boundingSphere = Sphere(center, glm::length(boundingCube.getBounds(1) - center));
	}

	VertexVis loadVertexVis(int v_id, int n_id, int t_id, tinyobj::attrib_t& attributes) {
		VertexVis vertexvis = {};
		vertexvis.pos = { attributes.vertices[3 * v_id + 0],
						  attributes.vertices[3 * v_id + 1],
						  attributes.vertices[3 * v_id + 2] };
		vertexvis.normal = { attributes.normals[3 * n_id + 0],
							 attributes.normals[3 * n_id + 1],
							 attributes.normals[3 * n_id + 2] };
		vertexvis.color = { 1,1,1 };
		vertexvis.tri_id = t_id;
		return vertexvis;
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

		for (tinyobj::shape_t& shape : shapes) primsize += shape.mesh.indices.size() / 3;

		vertices = std::vector<Vertex*>(attrib.vertices.size() / 3);
		indices = std::vector<glm::uint>(primsize * 3);
		verticesVis = std::vector<VertexVis>(primsize * 3);
		triangles = std::vector<Primitive*>(primsize);

		int h = 0;
		for (tinyobj::shape_t& shape : shapes) {
			for (int i = 0; i < shape.mesh.indices.size(); i += 3) {
				Primitive* prim = new Primitive({(h + i)/3});

				for (int j = 0; j < 3; j++) {
					tinyobj::index_t index = shape.mesh.indices[i + j];
					indices[h + i + j] = index.vertex_index;

					verticesVis[h + i + j] = loadVertexVis(index.vertex_index, index.normal_index, prim->id, attrib);
					
					if (vertices[index.vertex_index] == NULL) {
						Vertex* v = new Vertex({ index.vertex_index, verticesVis[h + i + j].pos });
						vertices[index.vertex_index] = v;
					}
					vertices[index.vertex_index]->triangles.push_back(prim);

					prim->vertices[j] = vertices[index.vertex_index];
				}

				prim->normal = glm::normalize(glm::cross(prim->vertices[0]->pos - prim->vertices[1]->pos, prim->vertices[2]->pos - prim->vertices[1]->pos));
				if (glm::dot(prim->normal, verticesVis[prim->id*3].normal) < 0) prim->normal = -prim->normal;
				prim->center = (prim->vertices[0]->pos + prim->vertices[1]->pos + prim->vertices[2]->pos) / 3.;

				for (int j = 0; j < 3; j++) {
					// Add the rays to the primitive

					int e_id = edges.size();
					Vertex* v0 = prim->vertices[j];
					Vertex* v1 = prim->vertices[(j + 1) % 3];
					Edge* e;
					if (v0->id < v1->id) e = new Edge{ e_id, { v0, v1 },  Ray(v0->pos, v1->pos, e_id) };
					else e = new Edge{ e_id, { v1, v0 },  Ray(v1->pos, v0->pos, e_id) };

					// Try inserting the edge to see if it already exists
					uint64_t edgeKey = (uint64_t)e->vertices[0]->id << 32 | (uint64_t)e->vertices[1]->id;
					if (!edges.contains(edgeKey)) {
						v0->edges.push_back(e);
						v1->edges.push_back(e);
						edges[edgeKey] = e;
					}
					else delete e;

					edges[edgeKey]->triangles.push_back(prim);
					prim->edges[j] = edges[edgeKey];
					prim->rays[j] = Ray(prim->vertices[j]->pos, prim->vertices[(j + 1) % 3]->pos, edges[edgeKey]->id);

				}

				triangles[prim->id] = prim;
			}
			h += shape.mesh.indices.size();
		}
		for (auto& e : edges) edgeVector.push_back(e.second);
	};

	RTCGeometry makeEmbreeGeom() {
		RTCGeometry model_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
		float* vertices_embr = (float*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), verticesVis.size());
		unsigned* triangles_embr = (unsigned*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), indices.size() / 3);
		for (int i = 0; i < verticesVis.size(); i++) {
			vertices_embr[i * 3] = verticesVis[i].pos.x;
			vertices_embr[i * 3 + 1] = verticesVis[i].pos.y;
			vertices_embr[i * 3 + 2] = verticesVis[i].pos.z;
		}
		for (int i = 0; i < indices.size(); i++)
			triangles_embr[i] = i;
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
		RTCGeometry model_geometry = makeEmbreeGeom();
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


	bool getIntersectionEmbree(const Ray& ray, int& primIndex, double& depth, bool backculling = false) {
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
				glm::vec3 norm = triangles[rayhit.hit.primID]->normal;
				if (glm::dot(norm, (glm::vec3)ray.direction) > 0) return false;
			}
			primIndex = rayhit.hit.primID;
			depth = rayhit.ray.tfar;
			return true;
		}
		//std::cout << "The ray did not collide with the object." << std::endl;

		return false;
	}

	//std::vector<Ray> getEdgeRaysForPrim(int prim) {
		//std::vector<Edge> triEdges = edgesPerTriangle[prim];
		//std::vector<Ray> edgeRays;
		//glm::vec3 normal = normalPerTri[prim];
		//glm::vec3 center = vertices[3 * prim].center;
		//Ray normalInvert = Ray(center + normal, center);
		//for (int i = 0; i < 3; i++) {
		//	Ray triEdge = Ray(verticesIndexed[triEdges[i].v[0]], verticesIndexed[triEdges[i].v[1]], triEdges[i].index);
		//	if (triEdge.side(normalInvert)) triEdge.inverseDir();
		//	edgeRays.push_back(triEdge);
		//}
		//return edgeRays;
	//}

	bool checkEdgeEdgeEdge(const Edge* e1, const Edge* e2, const Edge* e3) {
		uint64_t key;
		if (cacheEEE) {
			if (e1->id > e2->id) {
				if (e2->id > e3->id)		key = uint64_t(e1->id) << 42 | uint64_t(e2->id) << 21 | uint64_t(e3->id);
				else if (e1->id > e3->id)	key = uint64_t(e1->id) << 42 | uint64_t(e3->id) << 21 | uint64_t(e2->id);
				else						key = uint64_t(e3->id) << 42 | uint64_t(e1->id) << 21 | uint64_t(e2->id);
			}
			else {
				if (e3->id > e2->id)		key = uint64_t(e3->id) << 42 | uint64_t(e2->id) << 21 | uint64_t(e1->id);
				else if (e1->id > e3->id)	key = uint64_t(e2->id) << 42 | uint64_t(e1->id) << 21 | uint64_t(e3->id);
				else						key = uint64_t(e2->id) << 42 | uint64_t(e3->id) << 21 | uint64_t(e1->id);
			}
			if (edgeEdgeEdgeCombis.contains(key)) {
				cachehiteee++;
				return edgeEdgeEdgeCombis[key];
			}
		}

		std::vector<glm::dvec3> e1v = { e1->vertices[0]->pos, e1->vertices[1]->pos };
		std::vector<glm::dvec3> e2v = { e2->vertices[0]->pos, e2->vertices[1]->pos };
		Tetrahedron unboundTetra(e1v, e2v);

		std::vector<glm::dvec3> e3v = { e3->vertices[0]->pos, e3->vertices[1]->pos };
		if (unboundTetra.segmInTetra(e3v[0], e3v[1], 1E-7)) {
			if (cacheEEE) edgeEdgeEdgeCombis[key] = true;
			return true;
		}
		if (cacheEEE) edgeEdgeEdgeCombis[key] = false;
		return false;
	}

	bool checkEdgeEdgeCache(const Edge* e1, const Edge* e2, bool alldir, glm::vec3 maindir) {
		uint64_t key;
		int cacheIndex = 0;
		if (cacheEE) {
			if (!alldir) {
				int cacheIndex = maindir.x * 3 + maindir.y * 2 + maindir.z;
				cacheIndex < 0 ? cacheIndex += 6 : cacheIndex -= 3;
			}
			if (e1->id < e2->id) key = uint64_t(e1->id) << 32 | uint64_t(e2->id);
			else key = uint64_t(e2->id) << 32 | uint64_t(e1->id);
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

	bool checkEdgeEdge(const Edge* e1, const Edge* e2, bool alldir, glm::vec3 maindir) {
		//if (e1->vertices[0]->id == e2->vertices[0]->id || e1->vertices[1]->id == e2->vertices[1]->id ||
		//	e1->vertices[0]->id == e2->vertices[1]->id || e1->vertices[1]->id == e2->vertices[0]->id)  return false;

		glm::vec3 e1v1 = e1->vertices[0]->pos;
		glm::vec3 e1v2 = e1->vertices[1]->pos;
		glm::vec3 e2v1 = e2->vertices[0]->pos;
		glm::vec3 e2v2 = e2->vertices[1]->pos;
		if (!alldir) {
			if (!boundingCube.intersectSideBB(std::vector<Ray> {Ray(e1v1, e2v1), Ray(e1v1, e2v2), Ray(e1v2, e2v1), Ray(e1v2, e2v2)}, maindir)) return false;
		}

		// edge 1 to 2
		bool side1, side2;
		bool check12 = false;
		if (!e2->isSilhouetteForVertex(e1->vertices[0], side1)) {
			if (!e2->isSilhouetteForVertex(e1->vertices[1], side2)) 
				if (side1 == side2) return false;
		}

		// edge 2 to 1
		if (!e1->isSilhouetteForVertex(e2->vertices[0], side1)) {
			if (!e1->isSilhouetteForVertex(e2->vertices[1], side2)) 
				if (side1 == side2) return false;
		}
		return true;
	}

	bool checkSilhouetteEdge(Vertex * v, const Edge* e, bool alldir, glm::vec3& maindir, bool& side) {

		Vertex* v1 = e->vertices[0];
		Vertex* v2 = e->vertices[1];
		if (!alldir && (dot_fd(v1->pos, maindir) > dot_fd(v->pos, maindir) &&
			dot_fd(v2->pos, maindir) > dot_fd(v->pos, maindir))) return false;

		return e->isSilhouetteForVertex(v, side);

		//if (e->triangles.size() == 2) {
		//
		//	// define plane of triangle 1 and 2
		//	Plane p1(v1->pos, e->triangles[0]->normal);
		//	Plane p2(v2->pos, e->triangles[1]->normal);
		//	if (p1.pointOnPositiveSide(vpos) != p2.pointOnPositiveSide(vpos)) return 1;
		//	else { side = p1.pointOnPositiveSide(vpos); return 0; }
		//
		//	// optional check: if the triangle projected via vertices of edge falls outside box bounds
		//}
		//return true;
	}

	bool insertEdgeInMap(Edge* e, robin_hood::unordered_map<uint64_t, Edge*>& map) {
		uint64_t edgeKey = (uint64_t)e->vertices[0]->id << 32 | (uint64_t)e->vertices[1]->id;
		if (!map.contains(edgeKey)) {
			map[edgeKey] = e;
			return true;
		}
		return false;
	}

	//void findSilhouetteEdgesForTri(Primitive* tri, bool alldir, glm::vec3 maindir, std::vector<Edge*>& silhouetteEdge, std::set<int>& checktris = std::set<int>()) {

	//	robin_hood::unordered_map<uint64_t, Edge*> edgesToCheck;

	//	if (checktris.size() != 0) 	for (int t : checktris) for (Edge* e : triangles[t]->edges) insertEdgeInMap(e, edgesToCheck);
	//	else edgesToCheck = edges;

	//	glm::dvec3 invMaindir = glm::vec3(1.f) - maindir;
	//	glm::dvec3 maindirMin = boundingBox.getBounds(0);
	//	glm::dvec3 maindirMax = boundingBox.getBounds(1);
	//	Square window = boundingBox.getCubeSideSquare(maindir);

	//	for (auto& check_edge : edgesToCheck) {
	//		Edge* e = check_edge.second;
	//		// If edge is not in potential silhouettes
	//		if (!alldir && !e->silhouette) continue;

	//		// Don't use same triangle as silhouette edge.
	//		if (e->triangles[0] == tri) continue;
	//		if (e->triangles.size() == 2 && e->triangles[1] == tri) continue;

	//		// check if edge lies on correct side of triangle plane
	//		Plane triPlane = tri->getPlane();
	//		bool cntn = false;
	//		bool edgeIntersectPlane = true;
	//		for (Vertex* v : e->vertices) {
	//			if (!triPlane.pointOnPositiveSide(v->pos) || tri->hasVertex(v->id) || triPlane.pointOnPlane(v->pos, 1E-7)) {
	//				edgeIntersectPlane = true;
	//				cntn = true;
	//			}
	//			else cntn = false;
	//		}
	//		if (cntn) continue;

	//		for (Primitive* p : e->triangles) {
	//			if (glm::dot(window.normal, p->normal) < 0) {
	//				if (!window.intersectsPlaneFromLines(p->getRayVector())) cntn = true;
	//				else cntn = false;
	//			}
	//			else cntn = false;
	//		}
	//		if (cntn) continue;

	//		bool sideCheck = false;
	//		bool intersect = false;
	//		std::vector<glm::dvec3> intersectpts;
	//		bool first = true;
	//		bool silhouetteFound = false;

	//		glm::dvec3 v1pos = e->vertices[0]->pos;
	//		glm::dvec3 v2pos = e->vertices[1]->pos;

	//		for (int i = 0; i < 3; i++) {
	//			bool side = false;
	//			silhouetteFound = silhouetteFound || checkSilhouetteEdge(tri->vertices[i], e, alldir, maindir, side);
	//			if (!alldir) {
	//				Ray r1 = Ray(tri->vertices[i]->pos, v1pos);
	//				Ray r2 = Ray(tri->vertices[i]->pos, v2pos);

	//				intersect = intersect || window.inBounds(r1) || window.inBounds(r2);
	//				glm::vec3 windowR1Intersect = window.rayIntersection(r1);
	//				glm::vec3 windowR2Intersect = window.rayIntersection(r2);

	//				if (edgeIntersectPlane) {
	//					if (window.rayIntersectionDepth(r1) < 0) windowR1Intersect *= INFINITY;
	//					if (window.rayIntersectionDepth(r2) < 0) windowR2Intersect *= INFINITY;
	//				}
	//				intersectpts.push_back(windowR1Intersect);
	//				intersectpts.push_back(windowR2Intersect);
	//			}

	//			if (!silhouetteFound) {
	//				if (first) { first = false;  sideCheck = side; }
	//				else if (side != sideCheck) silhouetteFound = true;
	//			}

	//			if (silhouetteFound && (alldir || intersect)) { 
	//				silhouetteEdge.push_back(e); 
	//				break; 
	//			}
	//		}
	//		if (silhouetteFound && !intersect) {
	//			intersect = true;
	//			Cube bb(intersectpts);
	//			for (int i = 0; i < 3; i++) {
	//				if (invMaindir[i] == 1) {
	//					if (bb.getBounds(1)[i] < maindirMin[i] || bb.getBounds(0)[i] > maindirMax[i]) {
	//						intersect = false;
	//						break;
	//					}
	//				}
	//			}
	//			if (intersect) silhouetteEdge.push_back(e);
	//		}
	//	}
	//}



	bool getIntersectionNoAcceleration(const Ray& r, int& primIndex, float& depth) {
		primIndex = -1;
		depth = 10000.f;
		for (Primitive* prim : triangles) {
			double newdepth;
			if (prim->intersection(r, newdepth)) {
				if (newdepth < depth) {
					depth = newdepth;
					primIndex = prim->id;
				}
			}
		}
		return primIndex >= 0;
	};


	//void createIBO() {
	//	glGenVertexArrays(1, &vao2); // make VAO
	//	glBindVertexArray(vao2);
	//
	//	// copy interleaved vertex data (V/N/T) to VBO
	//   // GLuint vbo[2];
	//	GLuint vbo;
	//	glGenBuffers(1, &vbo);
	//	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	//	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW);
	//
	//	// copy index data to VBO
	//	GLuint ibo;
	//	glGenBuffers(1, &ibo);
	//	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	//	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);
	//
	//	// The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
	//	// Stride is the distance in bytes between vertices
	//	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	//	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
	//	glEnableVertexAttribArray(0);
	//
	//}

	void draw() {
		glBindVertexArray(vao);
		glDrawArrays(GL_TRIANGLES, 0, verticesVis.size());
	}

	void clearColors(glm::vec3& color) {
		//for (Vertex& v : vertices) v.color = color;
		for (int i = 0; i < verticesVis.size(); i++) verticesVis[i].color = color;
	}

	template <class T>
	void changePrimColors(T& setOfprims, glm::vec3& color) {
		for (auto i : setOfprims) changePrimColor(i, color);
	}

	void changePrimColor(int prim, glm::vec3& color) {
		for (int j = 0; j < 3; j++) {
			verticesVis[3 * prim + j].color = color;
		}
	}

	void changeSelected() {
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, verticesVis.size() * sizeof(VertexVis), verticesVis.data(), GL_STATIC_DRAW);
	}

	void createVAO() {
		//////////////////// Create Vertex Buffer Object
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, verticesVis.size() * sizeof(VertexVis), verticesVis.data(), GL_STATIC_DRAW);

		// Bind vertex data to shader inputs using their index (location)
		// These bindings are stored in the Vertex Array Object
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		// The position vectors should be retrieved from the specified Vertex Buffer Object with given offset and stride
		// Stride is the distance in bytes between vertices
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(VertexVis), reinterpret_cast<void*>(offsetof(VertexVis, pos)));
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(1, 3, GL_DOUBLE, GL_FALSE, sizeof(VertexVis), reinterpret_cast<void*>(offsetof(VertexVis, normal)));
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VertexVis), reinterpret_cast<void*>(offsetof(VertexVis, color)));
		glEnableVertexAttribArray(2);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(VertexVis), reinterpret_cast<void*>(offsetof(VertexVis, tri_id)));
		glEnableVertexAttribArray(3);
	}

};

#endif