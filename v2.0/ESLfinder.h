#pragma once
#include "cache.h"
#include "model.h"
#include "node.h"
#include "findLineThroughFour.h"
#include "combinations.h"
#include "rst.h"

class ESLfinder {
public:
	bool cache = false;
	Cache<std::vector<Line4>> combiCache;
	Model* model;
	bool alldir = false;
	glm::vec3 maindir;

	ESLfinder(bool cache) : cache(cache) {

	}

//	bool checkEdgeInLeafCombis(RaySpaceTree *rst, Edge& e, Node* leaf, std::vector<Ray>& splitLines, std::string combi_text, int nrOfsplitLines, std::vector<std::vector<int>>& combi, Ray& ray)
//	{
//		std::vector<std::vector<Ray>> edgeCombi;
//		int v1 = e.v[0];
//		int v2 = e.v[1];
//		glm::vec3 v1pos = model->verticesIndexed[v1];
//		glm::vec3 v2pos = model->verticesIndexed[v2];
//		Ray edge(model->verticesIndexed[v1], model->verticesIndexed[v2]);
//		if (nrOfsplitLines == 2) {
//			Ray v1ray = Ray(v1pos, v1pos + model->normalPerTri[e.triangles[0]]);
//			Ray v2ray = Ray(v2pos, v2pos + model->normalPerTri[e.triangles[0]]);
//			edgeCombi = { {edge, v1ray}, {edge, v2ray} };
//		}
//		else if (nrOfsplitLines == 3) edgeCombi = { {edge} };
//
//		for (std::vector<Ray>& rays : edgeCombi) {
//			for (int c = 0; c < combi.size(); c++) {
//				std::vector<Ray> lines;
//				for (int i = 0; i < nrOfsplitLines; i++) lines.push_back(splitLines[combi[c][i]]);
//				for (Ray& r : rays) lines.push_back(r);
//
//				std::vector<Line4> intersectLines = Lines4Finder::find(lines, model);
//				for (int i = 0; i < intersectLines.size() * 2; i++) {
//					ray = intersectLines[i / 2];
//					if (i % 2 == 1) ray.inverseDir();
//					if (!alldir && !checkRayInBox(ray, lines, nrOfsplitLines, false)) continue;
//					if (!checkRayInLeaf(rst, leaf, ray, lines, nrOfsplitLines, false)) continue;
//					if (nrOfsplitLines == 2) return true;
//					ray.get3DfromPlucker();
//					glm::vec3 r1, r2;
//					if (!model->boundingCube.intersect(ray, r1, r2)) continue;
//					glm::vec3 crossprod = glm::cross(v2pos - v1pos, r2 - r1);
//					float s = glm::dot(glm::cross(r1 - v1pos, r2 - r1), crossprod) / (glm::length(crossprod) * glm::length(crossprod));
//					float t = glm::dot(glm::cross(v1pos - r1, v2pos - v1pos), -crossprod) / (glm::length(crossprod) * glm::length(crossprod));
//
//					if (s >= 0.f && s <= 1.f && t >= 0.f && t <= 1.f) return true;
//				}
//			}
//		}
//		return false;
//	}
//
//	// add ignore list
//	bool checkRayInPrim(std::vector<Ray>& edgeRays, Ray& line, std::vector<Ray>& lines4, int prim, bool& inPlane, bool print) {
//		//check orientation
//		float orient = glm::dot((glm::vec3)line.direction, model->normalPerTri[prim]);
//		if (orient > 0) {
//			if (print) std::cout << "Ray not in Prim (wrong orientation): " << orient << std::endl << std::endl;
//			return false;
//		}
//
//		for (auto& er : edgeRays) {
//			bool equal = false;
//			for (auto& igray : lines4) {
//				if (igray.equal(er, 1E-8)) { //-8
//					equal = true;
//					break;
//				}
//			}
//			if (!equal && er.sideVal(line) < -1E-8) { //-8
//				if (print) std::cout << "Ray not in Prim (wrong side of edge): " << er.sideVal(line) << std::endl << std::endl;
//				return false;
//			}
//		}
//
//		if (line.inPlane(edgeRays[0], edgeRays[1], 1E-10)) {
//			if (print) std::cout << "Ray in plane: " << " ";
//			for (int i = 0; i < 3; i++) {
//				glm::dvec3 vert = model->vertices[3 * prim + i].pos;
//				if (line.throughVertex(vert, 1E-10)) {
//					if (print) std::cout << std::endl;
//					inPlane = true;
//					return true;
//				}
//			}
//			if (print) std::cout << std::endl;
//			return false;
//		}
//		return true;
//	}
//
//	bool checkRayInLeaf(RaySpaceTree* rst, Node* node, const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print) {
//
//		if (node == rst->rootNode) return true;
//		Node* parent = node->parent;
//		bool ign = false;
//		for (int i = 0; i < rayIgnoresize; i++) {
//			if (lines[i].equal(parent->splitter, 1E-8)) {
//				ign = true; //-8
//				break;
//			}
//		}
//		double sidev = parent->splitter.sideVal(ray);// (was other way around!) .sideVal;
//		if (abs(sidev) < 1E-8) ign = true;
//
//		if (parent->leftNode == node) {
//			if (sidev < 0 || ign) return checkRayInLeaf(rst, parent, ray, lines, rayIgnoresize, print);
//			if (print) std::cout << "Ray not in Leaf (on wrong side of splitting line, should be <0): " << sidev << std::endl << std::endl;
//		}
//		else {
//			if (sidev > 0 || ign) return checkRayInLeaf(rst, parent, ray, lines, rayIgnoresize, print);
//			if (print) std::cout << "Ray not in Leaf: on wrong side of splitting line, should be >0 " << sidev << std::endl << std::endl;
//		}
//		return false;
//	}
//
//	bool checkPrimVisibleForRay(Ray& ray, const int prim, std::vector<int>& ignore, std::vector<Ray>& lines, bool inplane, bool print) {
//
//		int embreePrim = -1;
//		float embreeDepth = 0.f;
//		float primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
//
//		if (ignore.size() > 0) {
//			std::vector<float> uniquePrimDepth;
//
//			// goes wrong right now
//			for (int i = 0; i < ignore.size() / 2; i++) {
//				float tryDepth = model->getIntersectionDepthForPrim(ignore[i * 2] * 3, ray);
//				bool found = false;
//				for (int j = 0; j < uniquePrimDepth.size(); j++) {
//					if (fabsf(tryDepth - uniquePrimDepth[j]) < 1E-3) {
//						found = true;
//						break;
//					}
//				}
//				if (!found) uniquePrimDepth.push_back(tryDepth);
//			}
//
//			int i = 0;
//			while (i < uniquePrimDepth.size()) {
//				//for (int i = 0; i < uniquePrimDepth.size(); i++) {
//				bool found = false;
//				model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
//				for (int j = 0; j < ignore.size() / 2; j++) {
//					// only checking first ignore primitive for now
//					float primdepth = model->getIntersectionDepthForPrim(ignore[j * 2] * 3, ray);
//					// check ray direction vs two triangles to see if it is still silhouette
//					if (ignore[j * 2 + 1] >= 0)
//						if ((glm::dot(model->normalPerTri[ignore[j * 2]], (glm::vec3)ray.direction) < 0) ==
//							(glm::dot(model->normalPerTri[ignore[j * 2 + 1]], (glm::vec3)ray.direction) < 0)) continue;
//
//					// Checking primary primdepth here does not result in 'correct' Extremal Stabbing Line
//					primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
//					//if (fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
//					if (fabsf(embreeDepth - primdepth) < 1E-3) {
//						glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
//						ray = Ray(neworig, neworig + ray.direction);
//						found = true;
//						break;
//					}
//				}
//				if (!found) {
//					if (lines.size() > 0) {
//						for (Ray& r : lines) {
//							glm::dvec3 intersectionTs = (ray.origin + ray.direction * (double)embreeDepth - r.origin) / r.direction;
//							if (fabsf(intersectionTs.x - intersectionTs.y) < 1E-10) {
//								glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
//								ray = Ray(neworig, neworig + ray.direction);
//								i--;
//								continue;
//							}
//						}
//					}
//					return false;
//				}
//				i++;
//			}
//		}
//
//		embreePrim = -1;
//		embreeDepth = 0.f;
//		primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
//		model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
//		if (embreePrim >= 0) {
//			if (embreePrim == prim || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
//			else if (inplane) {
//				glm::vec3 v1 = model->vertices[3 * embreePrim].pos;
//				glm::vec3 v2 = model->vertices[3 * embreePrim + 1].pos;
//				glm::vec3 v3 = model->vertices[3 * embreePrim + 2].pos;
//				std::vector<Ray> edgeRays = model->getEdgeRaysForPrim(prim);//{ Ray(v2, v1), Ray(v3, v2), Ray(v1, v3) };
//				glm::dvec3 cross = glm::cross(glm::normalize(edgeRays[0].direction), glm::normalize(edgeRays[1].direction));
//				double v = fabsf(glm::dot(ray.direction, cross));
//				if (v < 1E-10) return true;
//			}
//			else if (lines.size() > 0) {
//				for (Ray& r : lines) {
//					glm::dvec3 intersectionTs = (ray.origin + ray.direction * (double)embreeDepth - r.origin) / r.direction;
//					if (fabsf(intersectionTs.x - intersectionTs.y) < 1E-10) {
//						glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
//						ray = Ray(neworig, neworig + ray.direction);
//						primaryprimdepth = model->getIntersectionDepthForPrim(prim * 3, ray);
//						model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
//						if (embreePrim == prim || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
//					}
//				}
//			}
//		}
//		return false;
//	};
//
//	bool checkRayInBox(const Ray& ray, std::vector<Ray>& lines, int rayIgnoresize, bool print) {
//
//		if (!model->boundingCube.intersectSide(maindir, ray)) {
//			if (print) std::cout << "Ray not in Box (wrong side of edge)" << std::endl << std::endl;
//			return false;
//		}
//		return true;
//	}
//
//	bool checkCombi(RaySpaceTree* rst, const int prim, Line4& ray, Node* leaf, bool print, bool printAll, std::string combi_text, int combiNr,
//		int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges,
//		std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis,
//		std::vector<Ray>& silhouetteLines, std::vector<std::vector<int>>& silhLineCombis, std::vector<Edge>& silhouetteEdges,
//		std::vector<Ray>& silhVertexLines, std::vector<std::vector<int>>& silhVertexCombis,
//		std::vector<int>& silhouetteTris, std::vector<Ray>& triEdgeRays, bool vischeck = true) {
//
//		if (printAll) std::cout << combi_text << " combi's: " << combiNr << std::endl;
//
//		std::vector<std::vector<int>> edges = { {2,0}, {0,1}, {1,2} };
//		std::vector<int> vertexEdgeCheck;
//
//		for (int g = 0; g < std::max(1, std::min(nrOfTriEdges, 1) * 3); g++) {
//			if (nrOfTriEdges == 2) vertexEdgeCheck = { (int)model->indices[prim * 3 + edges[g][1]] };
//			if (nrOfTriEdges == 1) vertexEdgeCheck = { (int)model->indices[prim * 3 + edges[g][0]] , (int)model->indices[prim * 3 + edges[g][1]] };
//
//			for (int h = 0; h < std::max(1, (int)silhVertexCombis.size()); h++) {
//				if (nrOfVertices > 0) vertexEdgeCheck = { silhVertexCombis[h][2] };
//
//				for (int i = 0; i < std::max(1, (int)silhLineCombis.size()); i++) {
//					bool cntn = false;
//					for (int ix = 0; ix < nrOfsilhEdges; ix++) {
//						cntn = true;
//						if (nrOfsilhEdges > 0 && vertexEdgeCheck.size() == 1) {
//							Edge e = silhouetteEdges[silhLineCombis[i][ix]];
//							if (vertexEdgeCheck[0] == e.v[0] || vertexEdgeCheck[0] == e.v[1]) break;
//							bool side;
//							glm::vec3 vert = model->verticesIndexed[vertexEdgeCheck[0]];
//							if (model->checkSilhouetteEdge2(vert, e, true, glm::vec3(0), side) <= 0) break;
//							//extra check
//							//if (!model->boundingCube.intersectSideSwath(vert, model->verticesIndexed[e.v[0]], model->verticesIndexed[e.v[1]], maindir)) break;
//						}
//						//if ((nrOfsilhEdges > 2 || (nrOfsilhEdges == 2 && ix == 0)) && vertexEdgeCheck.size() > 0) {
//						//	std::vector<glm::vec3> n;
//						//	std::vector<float> d;
//						//	spaceSpannedByEdges(silhouetteEdges[silhLineCombis[i][ix]], silhouetteEdges[silhLineCombis[i][(ix + 1) % nrOfsilhEdges]], n, d);
//						//	std::vector<glm::vec3> checkPoints;
//						//	for (int vec : vertexEdgeCheck) checkPoints.push_back(model->verticesIndexed[vec]);
//						//	if (!checkPointsInHalfSpaces(n, d, checkPoints)) break;
//						//}
//						cntn = false;
//					}
//					if (cntn) continue;
//
//					for (int j = 0; j < std::max(1, (int)splitLineCombis.size()); j++) {
//
//						std::vector<Ray> lines;
//						std::vector<int> tris;
//						std::vector<uint64_t> indices;
//						int linecount = 0;
//						for (int k = 0; k < nrOfsplitLines; k++) {
//							lines.push_back(splitLines[splitLineCombis[j][k]]);
//							if (cache) indices.push_back(lines[linecount].index + model->edges.size() + model->vertices.size());
//							linecount++;
//						}
//						if (nrOfVertices) {
//							lines.push_back(silhouetteLines[silhVertexCombis[h][0]]);
//							lines.push_back(silhVertexLines[silhVertexCombis[h][1]]);
//							if (cache) {
//								indices.push_back(lines[linecount].index + model->edges.size());
//								indices.push_back(lines[linecount + 1].index + model->edges.size());
//							}
//							tris.push_back(silhouetteTris[2 * silhVertexCombis[h][0]]);
//							tris.push_back(silhouetteTris[2 * silhVertexCombis[h][0] + 1]);
//							linecount += 2;
//						}
//						for (int k = 0; k < nrOfsilhEdges; k++) {
//							lines.push_back(silhouetteLines[silhLineCombis[i][k]]);
//							tris.push_back(silhouetteTris[2 * silhLineCombis[i][k]]);
//							tris.push_back(silhouetteTris[2 * silhLineCombis[i][k] + 1]);
//							if (cache) indices.push_back(lines[linecount].index);
//							linecount++;
//						}
//						for (int k = 0; k < nrOfTriEdges; k++) {
//							lines.push_back(triEdgeRays[edges[g][k]]);
//							if (cache) indices.push_back(lines[linecount].index);
//							linecount++;
//						}
//
//						if (checkRaysThroughLines(rst, prim, ray, leaf, 0, triEdgeRays, printAll, tris, lines, indices, 4, vischeck)) { // change this!!
//							if (print) std::cout << combi_text << std::endl;
//							return true;
//						}
//					}
//				}
//			}
//		}
//		return false;
//	}
//
//
//	bool checkExtremalStabbingLine(RaySpaceTree* rst, const int prim, Ray& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
//		bool print, std::vector<int>& triIgnore, bool& inBox, std::vector<Ray>& lines,
//		int rayIgnoresize, bool vischeck) {
//
//		if (!alldir && ray.checkboth) {
//			if (!checkRayInBox(ray, lines, rayIgnoresize, print)) {
//				inBox = false;
//				return false;
//			}
//		}
//		if (!checkRayInLeaf(rst, leaf, ray, lines, rayIgnoresize, print)) return false;
//
//		ray.get3DfromPlucker();
//		bool inPlane = false;
//		if (!checkRayInPrim(edgeRays, ray, lines, prim, inPlane, print)) return false;
//
//		// not for generic viewing direction yet, only the three 'minima'
//		if (vischeck) {
//			double ttest = glm::dot(((model->boundingCube.getBounds(0) - glm::vec3(0.1) - (glm::vec3)ray.origin) / (glm::vec3)ray.direction), glm::abs(maindir));
//			glm::dvec3 neworig = glm::dvec3(ray.origin + ttest * ray.direction);
//			Ray ray2 = Ray(neworig, neworig + ray.direction);
//			if (!checkPrimVisibleForRay(ray2, prim, triIgnore, lines, inPlane, print)) return false;
//		}
//
//		return true;
//	}
//
//
//	bool checkRaysThroughLines(RaySpaceTree* rst, const int prim, Line4& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
//								bool print, std::vector<int>& visibleTriIgnore, std::vector<Ray>& lines,
//								std::vector<uint64_t>& indices, int rayIgnoresize, bool vischeck) {
//		std::vector<Line4> intersectLines;
//		uint64_t mapcombi;
//		if (!cache || !combiCache.getValue(indices, intersectLines)) intersectLines = Lines4Finder::find(lines, model);
//		std::vector<Line4> cacheLines;
//
//		for (int i = 0; i < intersectLines.size(); i++) {
//			ray = intersectLines[i];
//			bool inBox = true;
//			bool checkRay = checkRayAndReverse(rst, prim, ray, leaf, 0, edgeRays, print, visibleTriIgnore, inBox, lines, rayIgnoresize, vischeck);
//			if (cache && inBox) {
//				ray.checkboth = false;
//				cacheLines.push_back(ray);
//			}
//			if (checkRay) {
//				if (cache) {
//					if (i == 0) cacheLines.push_back(intersectLines[1]);
//					combiCache.storeValueAtLastKey(cacheLines);
//				}
//				return true;
//			}
//		}
//		if (cache) combiCache.storeValueAtLastKey(cacheLines);
//		return false;
//	}
//
//	bool checkRayAndReverse(RaySpaceTree* rst, const int prim, Line4& ray, Node* leaf, const int splitsize, std::vector<Ray>& edgeRays,
//		bool print, std::vector<int>& visibleTriIgnore, bool& inBox = false_bool, std::vector<Ray>& lines = std::vector<Ray>(),
//		int rayIgnoresize = 0, bool vischeck = true) {
//		if (checkExtremalStabbingLine(rst, prim, ray, leaf, 0, edgeRays, print, visibleTriIgnore, inBox, lines, rayIgnoresize, vischeck)) return true;
//		if (ray.checkboth && !inBox) {
//			ray.inverseDir();
//			inBox = true;
//			if (checkExtremalStabbingLine(rst, prim, ray, leaf, 0, edgeRays, print, visibleTriIgnore, inBox, lines, rayIgnoresize, vischeck)) return true;
//		}
//		return false;
//	}
//
//
//	bool checkVeVt(RaySpaceTree* rst, const int prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges,
//		std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//		if (printAll) std::cout << "V(e)V(t) combi's: " << silhouetteEdges.size() * 2 * 3 << std::endl;
//		// there are doubles in here! filter out to get unique vertices? a set perhaps??
//		for (int i = 0; i < silhouetteEdges.size(); i++) {
//			for (int v : silhouetteEdges[i].v) {
//				std::vector<int> tris = { silhouetteTris[i * 2], silhouetteTris[i * 2 + 1] };
//				for (int k = 0; k < 3; k++) {
//					ray = Line4(model->verticesIndexed[v], model->vertices[prim * 3 + k].pos);
//					if (checkRayAndReverse(rst, prim, ray, leaf, 0, edgeRays, printAll, tris)) {
//						if (print) std::cout << "V(e)V(t)" << std::endl;
//						return true;
//					}
//				}
//			}
//		}
//		return false;
//	}
//
//	bool checkVeVe(RaySpaceTree* rst, const int prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Edge>& silhouetteEdges,
//		std::vector<Ray>& silhouetteLines, std::vector<int>& silhouetteTris, std::vector<Ray>& edgeRays) {
//		//std::vector<std::vector<int>> combi = Combinations::combi2(silhouettesize);
//
//		if (printAll) std::cout << "V(e)V(e) combi's: " << silhouetteLines.size() * 2 * silhouetteLines.size() * 2 << std::endl;
//		std::vector<int> tris(4);
//		for (int i = 0; i < silhouetteEdges.size(); i++) {
//			Edge e1 = silhouetteEdges[i];
//			tris[0] = silhouetteTris[i * 2];
//			tris[1] = silhouetteTris[i * 2 + 1];
//			for (int v1 : e1.v) {
//				for (int j = i + 1; j < silhouetteEdges.size(); j++) {
//					Edge e2 = silhouetteEdges[j];
//					tris[2] = silhouetteTris[j * 2];
//					tris[3] = silhouetteTris[j * 2 + 1];
//
//					for (int v2 : e2.v) {
//						ray = Line4(model->verticesIndexed[v1], model->verticesIndexed[v2]);
//						if (checkRayAndReverse(rst, prim, ray, leaf, 0, edgeRays, printAll, tris)) {
//							if (print) std::cout << "V(e)V(e)" << std::endl;
//							//return true;
//						}
//					}
//				}
//			}
//		}
//		return false;
//	}
//
//	bool checkEdgeSplittingDuplicates(const int prim, std::vector<Ray>& splitLines, std::vector<Ray>& edgeRays, std::vector<bool>& sideLines) {
//		for (int i = 0; i < sideLines.size(); i++) {
//			for (int j = 0; j < edgeRays.size(); j++) {
//				if (splitLines[i].equal(edgeRays[j], 1E-5)) { // not very precise?!
//					glm::vec3 center = model->vertices[3 * prim].center;
//					glm::vec3 normal = model->normalPerTri[prim];
//					Ray checkRay = Ray(center + normal, center);
//
//					if (splitLines[i].side(checkRay) != sideLines[i]) {
//						return false;
//					}
//				}
//			}
//		}
//		return true;
//	}
//
//
//
//	bool checkSilhouetteCombis(RaySpaceTree* rst, const int prim, Line4& ray, Node* leaf, bool print, bool printAll, std::vector<Ray>& silhouetteLines,
//		std::vector<Ray>& splitLines, std::vector<Edge>& silhouetteEdgesToAdd, std::vector<Edge>& silhouetteEdges,
//		std::vector<int>& silhouetteTris, std::vector<Ray>& triEdgeRays,
//		std::set<int>& silhVertices, std::vector<Ray>& silhVertexRays, std::vector<std::vector<int>>& combiV) {
//
//
//		for (int i = 0; i < silhouetteEdgesToAdd.size(); i++) {
//			Edge toAdd = silhouetteEdgesToAdd[i];
//			silhouetteEdges.push_back(toAdd);
//			silhouetteLines.push_back(Ray(model->verticesIndexed[toAdd.v[0]], model->verticesIndexed[toAdd.v[1]], toAdd.index));
//			for (int t : toAdd.triangles) silhouetteTris.push_back(t);
//			if (toAdd.triangles.size() == 1) silhouetteTris.push_back(-1);
//			for (int v : toAdd.v) {
//				if (silhVertices.find(v) == silhVertices.end()) {
//					silhVertexRays.push_back(Ray(model->verticesIndexed[v], model->verticesIndexed[v] + model->normalPerTri[toAdd.triangles[0]], v));
//					combiV.push_back(std::vector<int>{ i, (int)silhVertices.size(), v });
//					silhVertices.insert(v);
//				}
//			}
//		}
//
//		std::vector<Ray> e1;
//		std::vector<std::vector<int>> e2;
//		std::vector<Edge> e3;
//		std::vector<int> e4;
//
//		if (checkVeVt(rst, prim, ray, leaf, print, printAll, silhouetteEdges, silhouetteTris, triEdgeRays)) return true;
//
//		std::vector<std::vector<int>> combi1S = Combinations::combi1(splitLines.size());
//
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SV(E)T", combi1S.size() * combiV.size() * 3, 1, 2, 0, 1,
//			splitLines, combi1S, silhouetteLines, e2, e3, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;
//
//		if (checkVeVe(rst, prim, ray, leaf, print, printAll, silhouetteEdges, silhouetteLines, silhouetteTris, triEdgeRays)) return true;
//
//		std::vector<std::vector<int>> combi1E = Combinations::combi1(silhouetteLines.size());
//
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SV(E)E", combi1S.size() * combi1E.size() * combiV.size(), 1, 2, 1, 0,
//			splitLines, combi1S, silhouetteLines, combi1E, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;
//
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SEV(T)", combi1S.size() * combi1E.size() * 3, 1, 0, 1, 2,
//			splitLines, combi1S, silhouetteLines, combi1E, silhouetteEdges, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//		std::vector<std::vector<int>> combi2S = Combinations::combi2(splitLines.size());
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSV(E)", combi2S.size() * combiV.size(), 2, 2, 0, 0,
//			splitLines, combi2S, silhouetteLines, e2, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;
//
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "V(E)ET", combi1E.size() * combiV.size() * 3, 0, 2, 1, 1,
//			e1, e2, silhouetteLines, combi1E, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;
//
//		std::vector<std::vector<int>> combi2E;
//		//getEECombis(rst, silhouetteEdges, silhouetteLines, combi2E);
//		if (combi2E.size() > 0) {
//			if (checkCombi(rst, prim, ray, leaf, print, printAll, "EEV(T)", combi2E.size() * 3, 0, 0, 2, 2,
//				e1, e2, silhouetteLines, combi2E, silhouetteEdges, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//			std::vector<std::vector<int>> combi3E;
//			//getEEECombis(rst, silhouetteEdges, silhouetteLines, combi2E, combi3E);
//
//			if (combi3E.size() > 0) {
//				if (checkCombi(rst, prim, ray, leaf, print, printAll, "EEET", combi3E.size() * 3, 0, 0, 3, 1,
//					e1, e2, silhouetteLines, combi3E, silhouetteEdges, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//				if (checkCombi(rst, prim, ray, leaf, print, printAll, "SEEE", combi1S.size() * combi3E.size(), 1, 0, 3, 0,
//					splitLines, combi1S, silhouetteLines, combi3E, e3, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//				std::vector<std::vector<int>> combi4E = Combinations::combiAddSelective(silhouetteLines.size(), combi3E);
//				if (combi4E.size() > 0) {
//					if (checkCombi(rst, prim, ray, leaf, print, printAll, "EEEE", combi4E.size() * 3, 0, 0, 4, 0,
//						e1, e2, silhouetteLines, combi4E, e3, e1, e2, silhouetteTris, triEdgeRays)) return true;
//				}
//			}
//
//			if (checkCombi(rst, prim, ray, leaf, print, printAll, "V(E)EE", combi2E.size() * combiV.size(), 0, 2, 2, 0,
//				e1, e2, silhouetteLines, combi2E, silhouetteEdges, silhVertexRays, combiV, silhouetteTris, triEdgeRays)) return true;
//
//			if (checkCombi(rst, prim, ray, leaf, print, printAll, "SEET", combi1S.size() * combi2E.size() * 3, 1, 0, 2, 1,
//				splitLines, combi1S, silhouetteLines, combi2E, silhouetteEdges, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//			if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSEE", combi2S.size() * combi2E.size(), 2, 0, 2, 0,
//				splitLines, combi2S, silhouetteLines, combi2E, e3, e1, e2, silhouetteTris, triEdgeRays)) return true;
//		}
//
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSET", combi2S.size() * combi1E.size() * 3, 2, 0, 1, 1,
//			splitLines, combi2S, silhouetteLines, combi1E, silhouetteEdges, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//		std::vector<std::vector<int>> combi3S = Combinations::combi3(splitLines.size());
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSSE", combi3S.size() * combi1E.size(), 3, 0, 1, 0,
//			splitLines, combi3S, silhouetteLines, combi1E, e3, e1, e2, silhouetteTris, triEdgeRays)) return true;
//
//		return false;
//	}
//
//
//	bool checkPrim(RaySpaceTree* rst, const int prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, std::vector<std::vector<int>>& combi4,
//		std::vector<Ray> splitLines, std::vector<bool>& sideLines, Line4& ray, Node* leaf, bool print, int edgeSelection, bool getedges,
//		std::vector<glm::vec3>& edges, std::vector<Ray>& eslEdges) {
//
//		std::vector<Ray> edgeRays = model->getEdgeRaysForPrim(prim);
//
//		bool printAll = false;
//		int splitLinesSize = splitLines.size();
//
//		// Check if triangle edges are same as splitting lines and if yes, if it lies on the correct side
//		if (!checkEdgeSplittingDuplicates(prim, splitLines, edgeRays, sideLines)) return false;
//
//		// Check basic combis of extremal stabbing lines
//		std::vector<Ray> e1;
//		std::vector<std::vector<int>> e2;
//		std::vector<Edge> e3;
//		std::vector<int> e4;
//
//		//if (!alldir && glm::dot(model->normalPerTri[prim], maindir) > 0) return false;
//
//		// Check without occlusion to see if prim can be excluded
//		if (!checkCombi(rst, prim, ray, leaf, false, printAll, "SSV(T)", combi2.size() * 3, 2, 0, 0, 2, splitLines, combi2, e1, e2, e3, e1, e2, e4, edgeRays, false) &&
//			!checkCombi(rst, prim, ray, leaf, false, printAll, "SSST", combi3.size() * 3, 3, 0, 0, 1, splitLines, combi3, e1, e2, e3, e1, e2, e4, edgeRays, false) &&
//			!checkCombi(rst, prim, ray, leaf, false, printAll, "SSSS", combi4.size(), 4, 0, 0, 0, splitLines, combi4, e1, e2, e3, e1, e2, e4, edgeRays, false)) return false;
//
//
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSV(T)", combi2.size() * 3, 2, 0, 0, 2, splitLines, combi2, e1, e2, e3, e1, e2, e4, edgeRays)) return true;
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSST", combi3.size() * 3, 3, 0, 0, 1, splitLines, combi3, e1, e2, e3, e1, e2, e4, edgeRays)) return true;
//
//		// Find silhouette edges for primitive
//		std::vector<Edge> silhouetteEdges;
//		std::vector<Edge> silhouetteEdgesFirst;
//		std::vector<Edge> silhouetteEdgesSecond;
//
//		if (edgeSelection == 0) {
//			model->findSilhouetteEdgesForTri(prim, alldir, maindir, silhouetteEdges);
//			//exact would be to test if triangles belong in leaf
//			for (int i = 0; i < silhouetteEdges.size(); i++) {
//				bool found = false;
//				for (int pr : silhouetteEdges[i].triangles) {
//					if (leaf->primitiveSet.find(pr) != leaf->primitiveSet.end()) {
//						silhouetteEdgesFirst.push_back(silhouetteEdges[i]);
//						found = true;
//						break;
//					}
//				}
//				if (!found) {
//					Ray eslEdge;
//					if (checkEdgeInLeafCombis(rst, silhouetteEdges[i], leaf, splitLines, "SSV(T)", 2, combi2, eslEdge) || //bug????
//						checkEdgeInLeafCombis(rst, silhouetteEdges[i], leaf, splitLines, "SSST", 3, combi3, eslEdge)) {
//						silhouetteEdgesSecond.push_back(silhouetteEdges[i]);
//						if (getedges) {
//							eslEdge.get3DfromPlucker();
//							eslEdges.push_back(eslEdge);
//						}
//					}
//				}
//			}
//		}
//		else if (edgeSelection == 1)
//			model->findSilhouetteEdgesForTri(prim, alldir, maindir, silhouetteEdgesFirst, leaf->parent->primitiveSet);
//
//		if (getedges) {
//			for (Edge e : silhouetteEdgesFirst) {
//				for (int v : e.v) edges.push_back(model->verticesIndexed[v]);
//			}
//			for (Edge e : silhouetteEdgesSecond) {
//				for (int v : e.v) edges.push_back(model->verticesIndexed[v]);
//			}
//		}
//
//		// Check all combis involving silhouette edges of some sort
//		int silhouettesize = 0;
//		std::vector<int> silhouetteTris;
//		std::vector<Ray> silhouetteLines;
//		std::set<int> silhEdgeVertices;
//		std::vector<Ray> silhVertexRays;
//		std::vector<std::vector<int>> edgeVertexCombis;
//		silhouetteEdges = std::vector<Edge>();
//
//		if (silhouetteEdgesFirst.size() > 0) {
//			if (checkSilhouetteCombis(rst, prim, ray, leaf, print, printAll, silhouetteLines, splitLines, silhouetteEdgesFirst, silhouetteEdges,
//				silhouetteTris, edgeRays, silhEdgeVertices, silhVertexRays, edgeVertexCombis)) return true;
//		}
//
//		// Check basic (but large) combi of extremal stabbing lines
//		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSSS", combi4.size(), 4, 0, 0, 0, splitLines, combi4, e1, e2, e3, e1, e2, e4, edgeRays)) return true;
//
//		// Check combis involving second tier silhouette edges
//		if (silhouetteEdgesSecond.size() > 0) {
//			if (checkSilhouetteCombis(rst, prim, ray, leaf, print, printAll, silhouetteLines, splitLines, silhouetteEdgesSecond, silhouetteEdges,
//				silhouetteTris, edgeRays, silhEdgeVertices, silhVertexRays, edgeVertexCombis)) return true;
//		}
//
//		return false;
//	}
//
};