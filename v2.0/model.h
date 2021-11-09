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
#include <set>
#include <utility>
#include "primitive.h"

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
	std::set<Edge, cmp_by_v> edges = std::set<Edge, cmp_by_v>();

	std::set<Edge, cmp_by_v> silhouetteEdges = std::set<Edge, cmp_by_v>();

	std::vector<std::vector<Edge>> edgesPerTriangle = std::vector<std::vector<Edge>>();
	std::vector<std::vector<Edge>> edgesPerVertex = std::vector<std::vector<Edge>>();

	std::vector<std::vector<int>> triPerVertex;
	std::vector<glm::vec3> normalPerTri = std::vector<glm::vec3>();

	int primsize = 0;
	Cube boundingBox;
	Cube boundingCube;
	glm::vec3 center = glm::vec3(0, 0, 0);
	float radius = 0.f;

	RTCDevice device;
	RTCScene scene;
	struct RTCIntersectContext context;

	Model() {};

	~Model() {
		rtcReleaseScene(scene);
	}

	Model(const char* filename, bool indexed = true) {
		loadModelFromFile(filename);
		findPotentialSilhouettes(glm::vec3(1, 0, 0));
		setUpEmbreeTracer();
		addGeometry(indexed);
		createVAO();
		createIBO();
		//enlargeModel();
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
		edgesPerVertex = std::vector<std::vector<Edge>>(attrib.vertices.size()/3);

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
						Edge e = Edge({ std::set<int>({verticespertri[x], verticespertri[(x + 1) % 3]}), std::vector<int>(), -1 });
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

		for (int i = 0; i < vertices.size(); i+=3) {
			glm::vec3 centroid = (vertices[i].pos + vertices[i+1].pos + vertices[i+2].pos) / 3.f;
			vertices[i].center = centroid;
			vertices[i + 1].center = centroid;
			vertices[i + 2].center = centroid;
		}

		boundingBox = Cube(glm::vec3(minx, miny, minz), glm::vec3(maxx, maxy, maxz));
		center = glm::vec3((minx + maxx) / 2.f, (miny + maxy) / 2.f, (minz + maxz) / 2.f);
		radius = glm::length(glm::vec3(maxx, maxy, maxz) - center);
		boundingCube = Cube(center, std::max(std::max(maxx-minx, maxy-miny),maxz-minz));
	};

	RTCGeometry makeEmbreeGeom(bool indexed = true) {
		RTCGeometry model_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
		int vertexsize = indexed ? vertices2.size() : vertices.size();
		float* vertices_embr = (float*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), vertexsize);
		unsigned* triangles_embr = (unsigned*)rtcSetNewGeometryBuffer(model_geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), indices.size() / 3);
		for (int i = 0; i < vertexsize; i++) {
			vertices_embr[i * 3] = indexed ? vertices2[i][0] : verticesL[i].pos.x;
			vertices_embr[i * 3 + 1] = indexed ? vertices2[i][1] : verticesL[i].pos.y;
			vertices_embr[i * 3 + 2] = indexed ? vertices2[i][2] : verticesL[i].pos.z;
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

	void findPotentialSilhouettes(glm::vec3 maindir) {
		std::vector<glm::vec3> cpoints = boundingCube.getCubeCornerPoints(maindir);
		for (auto &e : edges) {
			bool check = false;
			if (e.triangles.size() == 1) check = true;
			else {
				int count = 0;
				for (auto &c : cpoints) {
					glm::vec3 halfway = 0.5f * (vertices2[*e.vertices.begin()] + vertices2[*e.vertices.rbegin()]);
					bool dot1 = glm::dot(glm::normalize(halfway - c), normalPerTri[e.triangles[0]]) < 0;
					bool dot2 = glm::dot(glm::normalize(halfway - c), normalPerTri[e.triangles[1]]) < 0;
					if (dot1 != dot2) {
						check = true;
						break;
					}
					else if (!dot1) count++;
				}
				if (!check && count > 0 && count < 4) check = true;
			}
			if (check) {
				silhouetteEdges.insert(e);
				for (auto t : e.triangles) {
					vertices[3 * t].selected = -1.f;
					vertices[3 * t + 1].selected = -1.f;
					vertices[3 * t + 2].selected = -1.f;
				}
			}
		}
	}


	bool getIntersectionEmbree(const Ray& ray, int& primIndex, float& depth) {
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
			//glm::vec3 norm = vertices[rayhit.hit.primID * 3].normal;
			//if (glm::dot(norm, (glm::vec3)ray.direction) > 0) return false;
			primIndex = rayhit.hit.primID;
			depth = rayhit.ray.tfar;
			return true;
		}
		//std::cout << "The ray did not collide with the object." << std::endl;

		return false;
	}

	int checkSilhouetteEdge2(glm::vec3& vpos, const Edge& e, bool alldir, glm::vec3& maindir, bool& side) {

		if (e.triangles.size() == 2) {

			int v1 = *e.vertices.begin();
			int v2 = *e.vertices.rbegin();
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

			glm::vec3 v1 = vertices2[*e.vertices.begin()];
			glm::vec3 v2 = vertices2[*e.vertices.rbegin()];
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
			glm::vec3 jitter = glm::normalize(halfway - vpos- vertices[e.triangles[0] * 3].center);
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

		std::set<Edge, cmp_by_v> edgesToCheck = std::set<Edge, cmp_by_v>();

		if (checktris.size() != 0) 	for (int tri : checktris) for (Edge& edge : edgesPerTriangle[tri]) edgesToCheck.insert(edge);
		else edgesToCheck = edges;
		for (const Edge& e : edgesToCheck) {
			// If edge is not in potential silhouettes
			//if (!alldir && silhouetteEdges.find(e) == silhouetteEdges.end()) continue;

			// Don't use same triangle as silhouette edge.
			if (e.triangles[0] == prim || e.triangles[1] == prim) continue;

			bool sideCheck = false;
			for (int i = 0; i < 3; i++) {
				bool side = false;
				int checkEdge = checkSilhouetteEdge2(vertices[3 * prim + i].pos, e, alldir, maindir, side);
				if (checkEdge == 0) {
					if (i > 0 && side != sideCheck) silhouetteEdge.push_back(e);
					else { sideCheck = side; continue; }
				}
				else if (checkEdge > 0) silhouetteEdge.push_back(e);
				//break;
			}
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

	bool getIntersectionWithPrim(int i, const Ray& r, float &depth) {
		Ray edge;
		// clockwise triangle edges
		edge = Ray(vertices[i + 1].pos, vertices[i].pos);
		if (edge.side(r)) return false;
		edge = Ray(vertices[i + 2].pos, vertices[i + 1].pos);
		if (edge.side(r)) return false;
		edge = Ray(vertices[i].pos, vertices[i + 2].pos);
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
		glBufferData(GL_ARRAY_BUFFER, vertices2.size() * sizeof(glm::vec3), &vertices2[0], GL_STATIC_DRAW);

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
		glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, selected)));
		glEnableVertexAttribArray(4);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, id)));
		glEnableVertexAttribArray(5);
	}

};

#endif