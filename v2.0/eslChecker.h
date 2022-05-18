#pragma once
#include "rst.h"
#include "eslCandidate.h"

#define eps 1E-8

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
		
		if (esl.ray.checkboth) {
			if (!checkRayInWindow(esl.ray)) {
				esl.inBox = false;
				return false;
			}
		}
		if (!checkRayInLeaf(esl, leaf)) return false;

		if (primCheck && !checkRayInPrim(esl, visCheck)) return false;

		// not for generic viewing direction yet, only the three 'minima'
		if (visCheck && !checkPrimVisibleForRay(esl)) return false;

		return true;
	};

	bool checkRayInPrim(ESLCandidate& esl, bool visCheck) {
		if (esl.inPlane && esl.throughVertex) return true;
		//check orientation
		double orient = glm::dot(esl.ray.direction, prim->normal);
		if (orient > eps && !esl.inPlane) {
			if (print) std::cout << "Ray not in Prim (wrong orientation): " << orient << std::endl << std::endl;
			return false;
		}

		for (int i = 0; i < 3; i++) {
			if (!esl.doublesEdge(prim->edges[i])) {
				double side = prim->rays[i].sideVal(esl.ray);
				if (side < -eps) {
					if (print) std::cout << "Ray not in Prim (wrong side of edge): " << prim->rays[i].sideVal(esl.ray) << std::endl << std::endl;
					return false;
				}
				else if (esl.inPlane || (abs(orient) < 1E-8 && abs(side) < 1E-8)) {
					esl.inPlane = true;
					return checkRayInPlane(esl);
				}
			}
		}
		return true;
	};

	// add ray through triangle
	bool checkRayInPlane(ESLCandidate& esl) {
		if (print) std::cout << "Ray in plane " << " ";
		int leftSide = 0;
		if (esl.throughVertex) return true;
		for (Vertex* v : prim->vertices) {
			if (esl.ray.throughVertex(v, eps)) {
				if (print) std::cout << " and through vertex " << std::endl;
				return true;
			}
			if (Ray(v->pos, v->pos + prim->normal).side(esl.ray)) leftSide++;
		}
		if (leftSide > 0 && leftSide < 3) {
			esl.inPlaneNotVertex();
			return true;
		}
		if (print) std::cout << std::endl;
		return false;
	};

	bool checkRayInLeaf(ESLCandidate& esl, Node* node) {

		if (node == rst->rootNode) return true;
		Node* parent = node->parent;
		bool ign = esl.doublesSplit(parent->splitter);
		double sidev = parent->splitter.ray.sideVal(esl.ray);
		if (abs(sidev) < eps) ign = true;

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

	//bool checkSilhouettesForRay(ESLCandidate& esl, std::set<double>& intersectionDepths, double primaryprimdepth) {
	bool checkSilhouettesForRay(ESLCandidate & esl, std::set<int>& intersectionPrims, double primaryprimdepth) {

		for (Edge* e : esl.silhouetteEdges) {
			if (!e->isSilhouetteForRay(esl.ray)) return false;
			// check if ray intersects edge between vertices
			if (!e->intersectsRay(esl.ray)) return false;

			if (esl.ray.depthToIntersectionWithRay(e->ray) > primaryprimdepth) return false;
			// check if ray does not intersect exactly at vertex
			for (Vertex* v : e->vertices)
				if (esl.ray.throughPoint(v->pos, 1E-4)) return false;
			// store intersection depth
			//intersectionDepths.insert(esl.ray.depthToIntersectionWithRay(e->ray));
			for (Primitive* p : e->triangles) intersectionPrims.insert(p->id);//intersectionDepths.insert(p->getIntersectionDepth(esl.ray));
		}
		for (Vertex* v : esl.silhouetteVertices) {
			if (!esl.ray.throughPoint(v->pos, 1E-4)) return false;
			// check if line does not intersect mesh
			if (esl.ray.depthToPointOnRay(v->pos) > primaryprimdepth) return false;
			// store intersection depth
			//intersectionDepths.insert(esl.ray.depthToPointOnRay(v->pos));
			for (Primitive* p : v->triangles)  intersectionPrims.insert(p->id);// intersectionDepths.insert(p->getIntersectionDepth(esl.ray));

		}
		for (SplitSide& s : esl.splittingLines) {
			if (s.edge != NULL) { 
				glm::dvec3 pt = s.edge->ray.pointOfintersectWithRay(esl.ray);
				// check if intersection is on edge
				if (s.edge->pointOnRay(pt)) {
					// check if edge is silhouette
					double depthToPoint = esl.ray.depthToPointOnRay(pt);
					if (depthToPoint < primaryprimdepth) {
						if (s.edgeIsSilhouette) return false;
						// check if this edge is indeed a silhouetteEdge for this primitive (Edge* silhE : silhouette)
						if (!s.edge->isSilhouetteForRay(esl.ray)) return false;
						// if a ray interseting one of the silhouette triangles centers lies on the correct side of the
						// splitting line, the configuration is invalid.
						if (s.ray.side(Ray(s.edge->triangles[0]->center, s.edge->triangles[0]->center + esl.ray.direction)) == s.side) return false;

						for (Primitive* p : s.edge->triangles)  intersectionPrims.insert(p->id);// intersectionDepths.insert(p->getIntersectionDepth(esl.ray));
						//intersectionDepths.insert(depthToPoint);
					}
				}
			}
		}
		return true;
	};

	bool checkPrimVisibleForRay(ESLCandidate& esl) {
		
		double t;
		if (rst->alldir) {
			double tend;
			rst->model->boundingCube.intersect(esl.ray, t, tend);
		}
		else t = rst->window->rayIntersectionDepth(esl.ray) - 0.1f;
		esl.ray.offsetByDepth(t);

		std::set<int> intersectionPrims;
		double primaryprimdepth = prim->getIntersectionDepth(esl.ray, esl.inPlane, true);
		if (primaryprimdepth <= 0) primaryprimdepth = prim->getIntersectionDepth(esl.ray, true, true);

		if (!checkSilhouettesForRay(esl, intersectionPrims, primaryprimdepth)) return false;

		int embreePrim = -1;
		double embreeDepth = 0.f;

		// if not hit anything, return true (we already know it hits the prim)
		if (!rst->model->getIntersectionEmbree(esl.ray, embreePrim, embreeDepth)) return true;
		// if hit something, check it is same depth or further than prim
		double checkDepth = rst->model->triangles[embreePrim]->getIntersectionDepth(esl.ray);
		double depth = checkDepth;
		double offset = 0.001;
		double jump = offset + checkDepth;


		while (depth < primaryprimdepth - 1E-6 && embreePrim != prim->id) {
			bool inPlaneOfPrim = false;

			bool found = false;
			// check if it hits a silhouette vertex or edge
			for (int primid : intersectionPrims) {
				if (primid == embreePrim) {
					found = true;
					break;
				}
			}
			// check if primitive in same plane as target is hit
			if (esl.inPlane) {
				if (prim->getPlane().equal(rst->model->triangles[embreePrim]->getPlane())) {
					found = true;
					inPlaneOfPrim = true;
				}
			}

			if (!found) return false;

			// update the ray origin and shoot new ray
			esl.ray.offsetByDepth(jump);
			if (!rst->model->getIntersectionEmbree(esl.ray, embreePrim, embreeDepth)) return true;
			// update the depth with offset from previous intersection and new checkDepth
			Primitive* ePrim = rst->model->triangles[embreePrim];
			checkDepth = ePrim->getIntersectionDepth(esl.ray);
			if (fabs(checkDepth - embreeDepth) > offset || checkDepth < 0 || inPlaneOfPrim) {
				if (!ePrim->intersection(esl.ray, checkDepth, true)) jump = offset + embreeDepth;
				else {
					checkDepth = ePrim->getIntersectionDepth(esl.ray, true);
					jump = offset + checkDepth;
				}
			}
			else jump = offset + checkDepth;
			depth += jump;
		}
		return true;
	};

	bool checkRayInWindow(Line4& ray) {
		if (rst->alldir && !rst->model->boundingCube.intersect(ray)) {
			if (print) std::cout << "Ray not in Box" << std::endl << std::endl;
			return false;
		}
		else if (!rst->window->inBounds(ray, eps)) {
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