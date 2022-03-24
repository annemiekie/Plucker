#pragma once
#include "rst.h"
#include "eslCandidate.h"

class ESLChecker {
public:
	Node* leaf;
	Primitive* prim;
	RaySpaceTree* rst;

	bool print;

	ESLChecker() {};

	ESLChecker(Node* leaf, Primitive* prim, RaySpaceTree* rst, bool print = false) :
		leaf(leaf), prim(prim), rst(rst), print(print)
	{};

	bool isExtremalStabbing(ESLCandidate& esl, bool visCheck = true, bool primCheck = true) {
		return checkRayAndReverse(esl, visCheck, primCheck);
	};

	bool checkRayAndReverse(ESLCandidate& esl, bool visCheck, bool primCheck) {

		if (checkExtremalStabbingLine(esl, visCheck, primCheck)) return true;
		if ((rst->alldir && esl.inBox) || (!rst->alldir && esl.ray.checkboth && !esl.inBox)) {
			esl.ray.inverseDir();
			esl.inBox = true;
			if (checkExtremalStabbingLine(esl, visCheck, primCheck)) return true;
		}
		return false;
	};

	bool checkExtremalStabbingLine(ESLCandidate& esl, bool visCheck, bool primCheck) {
		
		esl.ray.get3DfromPlucker();
		if (esl.ray.checkboth) {
			if (!checkRayInBox(esl.ray)) {
				esl.inBox = false;
				return false;
			}
		}
		if (!checkRayInLeaf(esl, leaf)) return false;

		bool inPlane = false;
		if (primCheck && !checkRayInPrim(esl, inPlane)) return false;

		// not for generic viewing direction yet, only the three 'minima'
		if (visCheck && !checkPrimVisibleForRay(esl, inPlane)) return false;

		return true;
	};

	bool checkRayInPrim(ESLCandidate& esl, bool& inPlane) {
		//check orientation
		float orient = glm::dot(esl.ray.direction, prim->normal);
		if (orient > 0) {
			if (print) std::cout << "Ray not in Prim (wrong orientation): " << orient << std::endl << std::endl;
			return false;
		}
		for (auto& er : prim->rays) {
			bool equal = false;
			for (auto& igray : esl.lines4) {
				if (igray.equal(er, 1E-8)) {
					equal = true;
					break;
				}
			}
			if (!equal) {
				double side = er.sideVal(esl.ray);
				if (side < -1E-8) {
					if (print) std::cout << "Ray not in Prim (wrong side of edge): " << er.sideVal(esl.ray) << std::endl << std::endl;
					return false;
				}
				else if (abs(side) < 1E-10) return checkRayInPlane(esl.ray, inPlane);
			}
		}

		return true;
	};

	// add ray through triangle
	bool checkRayInPlane(Line4& ray, bool& inPlane) {
		if (print) std::cout << "Ray in plane " << " ";
		for (Vertex* v : prim->vertices) {
			if (ray.throughVertex(v, 1E-8)) {
				if (print) std::cout << " and through vertex " << std::endl;
				inPlane = true;
				return true;
			}
		}
		if (print) std::cout << std::endl;
		return false;
	};

	bool checkRayInLeaf(ESLCandidate& esl, Node* node) {

		if (node == rst->rootNode) return true;
		Node* parent = node->parent;
		bool ign = false;
		for (Ray& r : esl.lines4) {
			if (r.equal(parent->splitter.ray, 1E-8)) {
				ign = true;
				break;
			}
		}
		double sidev = parent->splitter.ray.sideVal(esl.ray);
		if (abs(sidev) < 1E-8) ign = true;

		if (parent->leftNode == node) {
			if (sidev < 0 || ign) return checkRayInLeaf(esl, parent);
			if (print) std::cout << "Ray not in Leaf (on wrong side of splitting line, should be <0): " << sidev << std::endl << std::endl;
		}
		else {
			if (sidev > 0 || ign) return checkRayInLeaf(esl, parent);
			if (print) std::cout << "Ray not in Leaf: on wrong side of splitting line, should be >0 " << sidev << std::endl << std::endl;
		}
		return false;
	};

	bool checkSilhouettesForRay(ESLCandidate& esl, std::set<double>& intersectionDepths, double primaryprimdepth) {
		for (Edge* e : esl.silhouetteEdges) {
			if (!e->isSilhouetteForRay(esl.ray)) return false;
			// check if ray intersects edge between vertices
			if (!e->intersectsRay(esl.ray)) return false;
			// check if ray does not intersect exactly at vertex
			for (Vertex* v : e->vertices)
				if (esl.ray.throughPoint(v->pos, 1E-4)) return false;
			// store intersection depth
			//intersectionDepths.insert(esl.ray.depthToIntersectionWithRay(e->ray));
			for (Primitive* p : e->triangles) intersectionDepths.insert(p->getIntersectionDepth(esl.ray));
		}
		for (Vertex* v : esl.silhouetteVertices) {
			if (!esl.ray.throughPoint(v->pos, 1E-4)) return false;
			// check if line does not intersect mesh
			// 
			// store intersection depth
			//intersectionDepths.insert(esl.ray.depthToPointOnRay(v->pos));
			for (Primitive* p : v->triangles) intersectionDepths.insert(p->getIntersectionDepth(esl.ray));

		}
		for (Split& s : esl.splittingLines) {
			if (s.edge != NULL) { 
				glm::dvec3 pt = s.edge->ray.pointOfintersectWithRay(esl.ray);
				// check if intersection is on edge
				if (s.edge->pointOnRay(pt)) {
					// check if edge is silhouette
					double depthToPoint = esl.ray.depthToPointOnRay(pt);
					if (depthToPoint < primaryprimdepth) {
						if (!s.edge->isSilhouetteForRay(esl.ray)) return false;
						else {
							for (Primitive* p : s.edge->triangles) intersectionDepths.insert(p->getIntersectionDepth(esl.ray));

							//intersectionDepths.insert(depthToPoint);
						}
					}
				}
			}
		}
		return true;
	};

	bool checkPrimVisibleForRay(ESLCandidate& esl, bool inplane) {
		
		double t;
		if (rst->alldir) {
			double tend;
			rst->model->boundingCube.intersect(esl.ray, t, tend);
		}
		else t = rst->model->boundingCube.getCubeSideSquare(rst->maindir).rayIntersectionDepth(esl.ray) - 0.1f;
		esl.ray.offsetByDepth(t);

		std::set<double> intersectionDepths;
		double primaryprimdepth = prim->getIntersectionDepth(esl.ray);
		if (!checkSilhouettesForRay(esl, intersectionDepths, primaryprimdepth)) return false;

		int embreePrim = -1;
		float embreeDepth = 0.f;

		// if not hit anything, return true (we already know it hits the prim)
		if (!rst->model->getIntersectionEmbree(esl.ray, embreePrim, embreeDepth)) return true;
		// if hit something, check it is same depth or further than prim
		double checkDepth = rst->model->triangles[embreePrim]->getIntersectionDepth(esl.ray);
		if (checkDepth > primaryprimdepth - 1E-6) return true;

		double depth = checkDepth;
		double offset = 0.001;

		while (depth < primaryprimdepth - 1E-6) {
			bool found = false;
			// check if it hits a silhouette vertex or edge
			for (double idepth : intersectionDepths) {
				if (fabsf(idepth - depth) < 1E-6) {
					found = true;
					intersectionDepths.erase(idepth);
					break;
				}
			}
			if (!found) return false;

			// update the ray origin and shoot new ray
			esl.ray.offsetByDepth(checkDepth + offset);
			if (!rst->model->getIntersectionEmbree(esl.ray, embreePrim, embreeDepth)) return true;
			// update the depth with offset from previous intersection and new checkDepth
			checkDepth = rst->model->triangles[embreePrim]->getIntersectionDepth(esl.ray);
			depth += (offset + checkDepth);
		}
		return true;
	};

	bool checkRayInBox(Line4& ray) {
		if (rst->alldir && !rst->model->boundingCube.intersect(ray)) {
			if (print) std::cout << "Ray not in Box" << std::endl << std::endl;
			return false;
		}
		else if (!rst->model->boundingCube.intersectSide(rst->maindir, ray)) {
			if (print) std::cout << "Ray not in Square" << std::endl << std::endl;
			return false;
		}
		return true;
	};

	// NEed to add to this
	bool checkNoSplitVertexProblem(Line4& ray) {
		// if ray through vertex
		bool throughvertex = false;
		for (Vertex* v : prim->vertices) if (ray.throughVertex(v)) throughvertex = true;
		if (!throughvertex) return true;

	}
};