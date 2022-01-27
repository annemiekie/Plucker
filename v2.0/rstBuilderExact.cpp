#include "rstBuilderExact.h"

void RSTBuilderExact::filterSplittingLines(RaySpaceTree *rst, Node* leaf, std::vector<Ray>& splitlines, std::vector<bool>& sides,
												std::vector<Ray>& filteredLines, std::vector<bool>& filteredSides) {
	std::vector<int> foundExtremalStabbing;
	std::vector<std::vector<int>>& splitCombi4 = Combinations::combi4(splitlines.size());
	//for (std::vector<int>& combi4 : splitCombi4) {
	for (int i = 0; i < splitCombi4.size(); i++) {
		std::vector<Ray> lines4;
		for (int c : splitCombi4[i]) lines4.push_back(splitlines[c]);
		std::vector<Line4> intersectLines = Lines4Finder::find(lines4);
		for (Line4& r : intersectLines) {
			if (checkRayInLeaf(rst, leaf, r, lines4)) {
				foundExtremalStabbing.push_back(i);
				break;
			}
			r.inverseDir();
			if (checkRayInLeaf(rst, leaf, r, lines4)) {
				foundExtremalStabbing.push_back(i);
				break;
			}
		}
	}
	std::set<int> filtered;
	for (int esl : foundExtremalStabbing) {
		for (int i : splitCombi4[esl]) {
			filtered.insert(i);
		}
	}
	for (int i : filtered) {
		filteredLines.push_back(splitlines[i]);
		if (i < sides.size()) filteredSides.push_back(sides[i]);
	}
}

// add ignore list
bool RSTBuilderExact::checkRayInPrim(Line4& ray, std::vector<Ray>& lines4, Primitive* prim, bool& inPlane, bool print) {
	//check orientation
	float orient = glm::dot((glm::vec3)ray.direction, prim->normal);
	if (orient > 0) {
		if (print) std::cout << "Ray not in Prim (wrong orientation): " << orient << std::endl << std::endl;
		return false;
	}

	for (auto& er : prim->rays) {
		bool equal = false;
		for (auto& igray : lines4) {
			if (igray.equal(er, 1E-8)) { //-8
				equal = true;
				break;
			}
		}
		if (!equal && er.sideVal(ray) < -1E-8) { //-8
			if (print) std::cout << "Ray not in Prim (wrong side of edge): " << er.sideVal(ray) << std::endl << std::endl;
			return false;
		}
	}

	if (ray.inPlane(prim->rays[0], prim->rays[1], 1E-10)) {
		if (print) std::cout << "Ray in plane: " << " ";
		for (int i = 0; i < 3; i++) {
			glm::dvec3 vert = prim->vertices[i]->pos;
			if (ray.throughVertex(vert, 1E-10)) {
				if (print) std::cout << std::endl;
				inPlane = true;
				return true;
			}
		}
		if (print) std::cout << std::endl;
		return false;
	}
	return true;
}

bool RSTBuilderExact::checkRayInLeaf(RaySpaceTree *rst, Node* node, const Line4& ray, std::vector<Ray>& lines4, bool print) {

	if (node == rst->rootNode) return true;
	Node* parent = node->parent;
	bool ign = false;
	for (Ray &r : lines4) {
		if (r.equal(parent->splitter, 1E-8)) {
			ign = true; //-8
			break;
		}
	}
	double sidev = parent->splitter.sideVal(ray);// (was other way around!) .sideVal;
	if (abs(sidev) < 1E-8) ign = true;

	if (parent->leftNode == node) {
		if (sidev < 0 || ign) return checkRayInLeaf(rst, parent, ray, lines4, print);
		if (print) std::cout << "Ray not in Leaf (on wrong side of splitting line, should be <0): " << sidev << std::endl << std::endl;
	}
	else {
		if (sidev > 0 || ign) return checkRayInLeaf(rst, parent, ray, lines4, print);
		if (print) std::cout << "Ray not in Leaf: on wrong side of splitting line, should be >0 " << sidev << std::endl << std::endl;
	}
	return false;
}

bool RSTBuilderExact::checkSilhouettesForRay(Line4& ray, std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges, 
												std::set<float>& intersectionDepths) {
	for (Edge* e : silhEdges) {
		if (!e->isSilhouetteForRay(ray)) return false;
		// check whether the edge is not intersected at the vertex
		for (Vertex* v : e->vertices)
			if (ray.throughPoint(v->pos, 1E-8)) return false;
		// store intersection depth
		intersectionDepths.insert(ray.depthToIntersectionWithRay(e->ray));
	}
	for (Vertex* v : silhVertices) {
		// check if line goes through vertex
		if (!ray.throughPoint(v->pos, 1E-8)) return false;
		// check if line does not intersect mesh

		// store intersection depth
		intersectionDepths.insert(ray.depthToPointOnRay(v->pos));
	}
}

bool RSTBuilderExact::checkPrimVisibleForRay(RaySpaceTree* rst, Line4& ray, Primitive* prim, 
											std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges,
											std::vector<Ray>& lines4, bool inplane, bool print) {

	float t;
	if (rst->alldir) {
		float tend;
		rst->model->boundingCube.intersect(ray, t, tend);
	}
	else t = rst->model->boundingCube.getCubeSideSquare(rst->maindir).rayIntersectionDepth(ray) - 0.1f;
	glm::dvec3 neworig = glm::dvec3(ray.origin + (double)t * ray.direction);
	ray = Line4(neworig, neworig + ray.direction);

	std::set<float> intersectionDepths;
	if (!checkSilhouettesForRay(ray, silhVertices, silhEdges, intersectionDepths)) return false;

	int embreePrim = -1;
	float embreeDepth = 0.f;

	std::set<float>::iterator intDepth = intersectionDepths.begin();
	while (intDepth != intersectionDepths.end()) {
		rst->model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
		bool found = false;
		if (fabsf(embreeDepth - *intDepth) < 1E-3) { found = true; intDepth++; }
		else {
			for (Ray& r : lines4)
				if (ray.intersectsWithRayAtDepth(r, embreeDepth)) { found = true; break; }
		}
		if (found) {
			glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
			ray = Line4(neworig, neworig + ray.direction);
		}
		else return false;
	}

	embreePrim = -1;
	embreeDepth = 0.f;
	float primaryprimdepth = prim->getIntersectionDepth(ray);
	rst->model->getIntersectionEmbree(ray, embreePrim, embreeDepth);

	if (embreePrim >= 0) {
		if (embreePrim == prim->id || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
		else if (inplane) return true;
		else if (lines4.size() > 0) {
			for (Ray& r : lines4) {
				if (ray.intersectsWithRayAtDepth(r, embreeDepth)) {
					glm::dvec3 neworig = glm::dvec3(ray.origin + (embreeDepth + 0.001) * ray.direction);
					ray = Line4(neworig, neworig + ray.direction);
					primaryprimdepth = prim->getIntersectionDepth(ray);
					rst->model->getIntersectionEmbree(ray, embreePrim, embreeDepth);
					if (embreePrim == prim->id || fabsf(embreeDepth - primaryprimdepth) < 1E-3) return true;
				}
			}
		}
	}
	return false;
};

bool RSTBuilderExact::checkRayInBox(RaySpaceTree* rst, const Line4& ray, std::vector<Ray>& lines, bool print) {
	if (rst->alldir && !rst->model->boundingCube.intersect(ray)) {
		if (print) std::cout << "Ray not in Box" << std::endl << std::endl;
		return false;
	}
	else if (!rst->model->boundingCube.intersectSide(rst->maindir, ray)) {
		if (print) std::cout << "Ray not in Square" << std::endl << std::endl;
		return false;
	}
	return true;
}

bool RSTBuilderExact::checkCombi(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, 
									std::string combi_text, int combiNr, 
									int nrOfsplitLines, int nrOfVertices, int nrOfsilhEdges, int nrOfTriEdges, 
									std::vector<Ray>& splitLines, std::vector<std::vector<int>>& splitLineCombis, bool vischeck,
									std::vector<std::vector<int>>& silhLineCombis, std::vector<Edge*>& silhouetteEdges, 
									std::vector<Vertex*> silhouetteVertices) {

	if (printAll) std::cout << combi_text << " combi's: " << combiNr << std::endl;

	std::vector<Vertex*> vertexEdgeCheck;

	for (int g = 0; g < std::max(1, std::min(nrOfTriEdges, 1) * 3); g++) {
		if (nrOfTriEdges == 2) vertexEdgeCheck = { prim->vertices[g] };
		//if (nrOfTriEdges == 1) vertexEdgeCheck = { (int)model->indices[prim * 3 + edges[g][0]] , (int)model->indices[prim * 3 + edges[g][1]] };

		for (int h = 0; h < std::max(1, (int)silhouetteVertices.size()); h++) {
			if (nrOfVertices > 0) vertexEdgeCheck = { silhouetteVertices[h] };

			for (int i = 0; i < std::max(1, (int)silhLineCombis.size()); i++) {
				bool cntn = false;
				for (int ix = 0; ix < nrOfsilhEdges; ix++) {
					cntn = true;
					if (nrOfsilhEdges > 0 && vertexEdgeCheck.size() == 1) {
						Edge* e = silhouetteEdges[silhLineCombis[i][ix]];
						// if silhouette edge vertices are the same as one of the vertices used break
						if (vertexEdgeCheck[0] == e->vertices[0] || vertexEdgeCheck[0] == e->vertices[1]) break;
						//glm::vec3 vert = rst->model->verticesIndexed[vertexEdgeCheck[0]];
						bool side;
						if (!e->isSilhouetteForPos(vertexEdgeCheck[0]->pos, side)) break;
						// rst->model->checkSilhouetteEdge(vert, e, true, glm::vec3(0), side) <= 0) break;
						//extra check
						//if (!model->boundingCube.intersectSideSwath(vert, model->verticesIndexed[e.v[0]], model->verticesIndexed[e.v[1]], maindir)) break;
					}
					//if ((nrOfsilhEdges > 2 || (nrOfsilhEdges == 2 && ix == 0)) && vertexEdgeCheck.size() > 0) {
					//	std::vector<glm::vec3> n;
					//	std::vector<float> d;
					//	spaceSpannedByEdges(silhouetteEdges[silhLineCombis[i][ix]], silhouetteEdges[silhLineCombis[i][(ix + 1) % nrOfsilhEdges]], n, d);
					//	std::vector<glm::vec3> checkPoints;
					//	for (int vec : vertexEdgeCheck) checkPoints.push_back(model->verticesIndexed[vec]);
					//	if (!checkPointsInHalfSpaces(n, d, checkPoints)) break;
					//}
					cntn = false;
				}
				if (cntn) continue;

				for (int j = 0; j < std::max(1, (int)splitLineCombis.size()); j++) {

					std::vector<Ray> lines4;
					std::vector<Edge*> silhEdges;
					std::vector<Vertex*> silhVertices;
					std::vector<uint64_t> indices;
					int linecount = 0;
					for (int k = 0; k < nrOfsplitLines; k++) {
						lines4.push_back(splitLines[splitLineCombis[j][k]]);
					//	if (cacheCombi) indices.push_back(lines[linecount].index + model->edges.size() + model->vertices.size());
						linecount++;
					}
					if (nrOfVertices) {
						lines4.push_back(silhouetteEdges[h / 2]->vertices[h % 2]->edges[0]->ray);
						lines4.push_back(silhouetteEdges[h / 2]->vertices[h % 2]->edges[1]->ray);
						silhVertices = { silhouetteEdges[h / 2]->vertices[h % 2] };
					//	if (cacheCombi) {
					//		indices.push_back(lines[linecount].index + model->edges.size());
					//		indices.push_back(lines[linecount + 1].index + model->edges.size());
					//	}
						linecount += 2;
					}
					for (int k = 0; k < nrOfsilhEdges; k++) {
						lines4.push_back(silhouetteEdges[silhLineCombis[i][k]]->ray);
						silhEdges.push_back(silhouetteEdges[silhLineCombis[i][k]]);
					//	if (cacheCombi) indices.push_back(lines[linecount].index);
						linecount++;
					}
					for (int k = 0; k < nrOfTriEdges; k++) {
						lines4.push_back(prim->rays[(g+k)%3]);// triEdgeRays[edges[g][k]]);
					//	if (cacheCombi) indices.push_back(lines[linecount].index);
						linecount++;
					}

					if (checkRaysThroughLines(rst, prim, ray, leaf, printAll, lines4, indices, vischeck)) { // change this!!
						if (print) std::cout << combi_text << std::endl;
						return true;
					}
				}
			}
		}
	}
	return false;
}


bool RSTBuilderExact::checkExtremalStabbingLine(RaySpaceTree *rst, Primitive* prim, Line4& ray, Node* leaf, 
												bool print, bool& inBox, std::vector<Ray>& lines,
												bool vischeck, std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges) {

	if (ray.checkboth) {
		if (!checkRayInBox(rst, ray, lines, print)) {
			inBox = false;
			return false;
		}
	}
	if (!checkRayInLeaf(rst, leaf, ray, lines, print)) return false;

	ray.get3DfromPlucker();
	bool inPlane = false;
	if (!checkRayInPrim(ray, lines, prim, inPlane, print)) return false;

	// not for generic viewing direction yet, only the three 'minima'
	if (vischeck && !checkPrimVisibleForRay(rst, ray, prim, silhVertices, silhEdges, lines, inPlane, print)) return false;

	return true;
}

bool RSTBuilderExact::checkRaysThroughLines(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, 
											std::vector<Ray>& lines, std::vector<uint64_t>& indices, bool vischeck,
											std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges) {

	std::vector<Line4> intersectLines;
	//uint64_t mapcombi;
	//if (!cacheCombi || !combiCache.getValue(indices, intersectLines)) intersectLines = LineThroughFour::find(lines, model);
	std::vector<Ray> cacheLines;

	for (int i = 0; i < intersectLines.size(); i++) {
		ray = intersectLines[i];
		bool inBox = true;
		bool checkRay = checkRayAndReverse(rst, prim, ray, leaf, print, inBox, lines, vischeck, silhVertices, silhEdges);
		//if (cacheCombi && inBox) {
		//	ray.checkboth = false;
		//	cacheLines.push_back(ray);
		//}
		if (checkRay) {
			//if (cacheCombi) {
			//	if (i == 0) cacheLines.push_back(intersectLines[1]);
			//	combiCache.storeValueAtLastKey(cacheLines);
			//}
			return true;
		}
	}
	//if (cacheCombi) combiCache.storeValueAtLastKey(cacheLines);
	return false;
}

bool RSTBuilderExact::checkRayAndReverse(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, 										
										bool& inBox, std::vector<Ray>& lines, bool vischeck,
										std::vector<Vertex*>& silhVertices, std::vector<Edge*>& silhEdges) {

	if (checkExtremalStabbingLine(rst, prim, ray, leaf, print, inBox, lines, vischeck, silhVertices, silhEdges)) return true;
	if ((rst->alldir && inBox) || (!rst->alldir && ray.checkboth && !inBox)) {
		ray.inverseDir();
		inBox = true;
		if (checkExtremalStabbingLine(rst, prim, ray, leaf, print, inBox, lines, vischeck, silhVertices, silhEdges)) return true;
	}
	return false;
}


bool RSTBuilderExact::checkVeVt(RaySpaceTree *rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, 
								std::vector<Vertex*>& silhVertices) {
	if (printAll) std::cout << "V(e)V(t) combi's: " << silhVertices.size() * 3 << std::endl;
	for (Vertex* vs : silhVertices) {
		for (Vertex* vt : prim->vertices) {
			ray = Line4(vs->pos, vt->pos);
			bool x;
			std::vector<Vertex*> vertices = { vs };
			if (checkRayAndReverse(rst, prim, ray, leaf, printAll, x, std::vector<Ray>(), true, vertices)) {
				if (print) std::cout << "V(e)V(t)" << std::endl;
				return true;
			}
		}
	}
	return false;
}

bool RSTBuilderExact::checkVeVe(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, bool printAll, 
								std::vector<Vertex*>& silhVertices) {
	//std::vector<std::vector<int>> combi = Combinations::combi2(silhouettesize);

	if (printAll) std::cout << "V(e)V(e) combi's: " << silhVertices.size() *(silhVertices.size()-1) << std::endl;
	for (Vertex* vs : silhVertices) {
		for (Vertex* vt : silhVertices) {
			if (vs->id == vt->id) continue;
			ray = Line4(vs->pos, vt->pos);
			bool x;
			std::vector<Vertex*> vertices = { vs, vt };
			if (checkRayAndReverse(rst, prim, ray, leaf, printAll, x, std::vector<Ray>(), true, vertices)) {
				if (print) std::cout << "V(e)V(e)" << std::endl;
				return true;
			}
		}
	}
	return false;
}

bool RSTBuilderExact::checkEdgeInLeafCombis(RaySpaceTree *rst, Edge* e, Node* leaf, std::vector<Ray>& splitLines, 
												std::vector<std::vector<int>>& combi, Line4& ray)
{
	std::vector<std::vector<Ray>> edgeCombi;
	if (combi[0].size() == 2) edgeCombi = { {e->vertices[0]->edges[0]->ray, e->vertices[0]->edges[1]->ray},
											{e->vertices[1]->edges[0]->ray, e->vertices[1]->edges[1]->ray} };
	else if (combi[0].size() == 3) edgeCombi = { {e->ray} };
	else return false;

	for (std::vector<Ray>& rays : edgeCombi) {
		for (int c = 0; c < combi.size(); c++) {
			std::vector<Ray> lines;
			for (int i = 0; i < combi[0].size(); i++) lines.push_back(splitLines[combi[c][i]]);
			for (Ray& r : rays) lines.push_back(r);

			std::vector<Line4> intersectLines = Lines4Finder::find(lines);

			for (int i = 0; i < intersectLines.size() * 2; i++) {
				ray = intersectLines[i / 2];
				if (i % 2 == 1) ray.inverseDir();
				if (!checkRayInBox(rst, ray, lines, false)) continue;
				if (!checkRayInLeaf(rst, leaf, ray, lines, false)) continue;
				if (combi[0].size() == 2) return true;

				ray.get3DfromPlucker();
				float t = e->ray.depthToIntersectionWithRay(ray);
				if (t > 0 && t < 1) return true;

				/*glm::vec3 r1, r2;
				glm::vec3 crossprod = glm::cross(v2pos - v1pos, r2 - r1);
				float s = glm::dot(glm::cross(r1 - v1pos, r2 - r1), crossprod) / (glm::length(crossprod) * glm::length(crossprod));
				float t = glm::dot(glm::cross(v1pos - r1, v2pos - v1pos), -crossprod) / (glm::length(crossprod) * glm::length(crossprod));

				if (s >= 0.f && s <= 1.f && t >= 0.f && t <= 1.f) return true;*/
			}
		}
	}
	return false;
}


void RSTBuilderExact::getEEECombis(RaySpaceTree* rst, std::vector<Edge*>& silhouetteEdges,
									std::vector<std::vector<int>>& combi2Edges, std::vector<std::vector<int>>& combi3Edges) {
	if (silhouetteEdges.size() >= 3 && combi2Edges.size() > 0) {
		std::vector<std::vector<int>> combi = Combinations::combiAddSelective(silhouetteEdges.size(), combi2Edges);
		for (std::vector<int>& c : combi) {
			Edge* e1 = silhouetteEdges[c[0]];
			Edge* e2 = silhouetteEdges[c[1]];
			Edge* e3 = silhouetteEdges[c[2]];

			if (rst->model->checkEdgeEdgeEdge(e1, e2, e3)) combi3Edges.push_back(c);
		}
	}
}

void RSTBuilderExact::getEECombis(RaySpaceTree* rst, std::vector<Edge*>& silhouetteEdges, std::vector<std::vector<int>>& combi2Edges) {
	if (silhouetteEdges.size() >= 2) {
		std::vector<std::vector<int>> combi = Combinations::combi2(silhouetteEdges.size());

		for (std::vector<int>& c : combi) {
			Edge* e1 = silhouetteEdges[c[0]];
			Edge* e2 = silhouetteEdges[c[1]];
			if (rst->model->checkEdgeEdgeCache(e1, e2, rst->alldir, rst->maindir)) combi2Edges.push_back(c);
		}
	}
}


bool RSTBuilderExact::checkSilhouetteCombis(RaySpaceTree *rst, Primitive *prim, Line4& ray, Node* leaf, bool print, bool printAll,
										std::vector<Ray>& splitLines, std::vector<Edge*>& silhouetteEdgesToAdd, std::vector<Edge*>& silhouetteEdges,
										std::set<Vertex*, Vertex::cmp_ptr>& silhVertices) {


for (Edge* toAdd : silhouetteEdgesToAdd) {
	silhouetteEdges.push_back(toAdd);
	for (Vertex* v : toAdd->vertices) silhVertices.insert(v);
	//{
	//	if (silhVertices.find(v) == silhVertices.end()) {
	//		silhVertexRays.push_back(Ray(model->verticesIndexed[v], model->verticesIndexed[v] + model->normalPerTri[toAdd.triangles[0]], v));
	//		combiV.push_back(std::vector<int>{ i, (int)silhVertices.size(), v });
	//		silhVertices.insert(v);
	//	}
	//}
}

std::vector<Vertex*> silhouetteVertices;
for (Vertex* v : silhVertices) silhouetteVertices.push_back(v);

std::vector<Ray> e1;
std::vector<std::vector<int>> e2;
std::vector<Edge*> e3;
std::vector<int> e4;

if (checkVeVt(rst, prim, ray, leaf, print, printAll, silhouetteVertices)) return true;

std::vector<std::vector<int>> combi1S = Combinations::combi1(splitLines.size());

if (checkCombi(rst, prim, ray, leaf, print, printAll, "SV(E)T", combi1S.size() * silhouetteEdges.size() * 2 * 3, 1, 2, 0, 1,
				splitLines, combi1S, true, e2, e3, silhouetteVertices)) return true;

if (checkVeVe(rst, prim, ray, leaf, print, printAll, silhouetteVertices)) return true;

std::vector<std::vector<int>> combi1E = Combinations::combi1(silhouetteEdges.size());

if (checkCombi(rst, prim, ray, leaf, print, printAll, "SV(E)E", combi1S.size() * combi1E.size() * silhouetteEdges.size() * 2, 1, 2, 1, 0,
				splitLines, combi1S, true, combi1E, silhouetteEdges, silhouetteVertices)) return true;

if (checkCombi(rst, prim, ray, leaf, print, printAll, "SEV(T)", combi1S.size() * combi1E.size() * 3, 1, 0, 1, 2,
				splitLines, combi1S, true, combi1E, silhouetteEdges)) return true;

std::vector<std::vector<int>> combi2S = Combinations::combi2(splitLines.size());
if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSV(E)", combi2S.size() * silhouetteEdges.size() * 2, 2, 2, 0, 0,
				splitLines, combi2S, true, e2, e3, silhouetteVertices)) return true;

if (checkCombi(rst, prim, ray, leaf, print, printAll, "V(E)ET", combi1E.size() * silhouetteEdges.size() * 2 * 3, 0, 2, 1, 1,
				e1, e2, true, combi1E, silhouetteEdges, silhouetteVertices)) return true;

std::vector<std::vector<int>> combi2E;
getEECombis(rst, silhouetteEdges, combi2E);
if (combi2E.size() > 0) {
	if (checkCombi(rst, prim, ray, leaf, print, printAll, "EEV(T)", combi2E.size() * 3, 0, 0, 2, 2,
					e1, e2, true, combi2E, silhouetteEdges)) return true;

	std::vector<std::vector<int>> combi3E;
	getEEECombis(rst, silhouetteEdges, combi2E, combi3E);

	if (combi3E.size() > 0) {
		if (checkCombi(rst, prim, ray, leaf, print, printAll, "EEET", combi3E.size() * 3, 0, 0, 3, 1,
						e1, e2, true, combi3E, silhouetteEdges)) return true;

		if (checkCombi(rst, prim, ray, leaf, print, printAll, "SEEE", combi1S.size() * combi3E.size(), 1, 0, 3, 0,
						splitLines, combi1S, true, combi3E, silhouetteEdges)) return true;

		std::vector<std::vector<int>> combi4E = Combinations::combiAddSelective(silhouetteEdges.size(), combi3E);
		if (combi4E.size() > 0) {
			if (checkCombi(rst, prim, ray, leaf, print, printAll, "EEEE", combi4E.size() * 3, 0, 0, 4, 0,
							e1, e2, true, combi4E, silhouetteEdges)) return true;
		}
	}

	if (checkCombi(rst, prim, ray, leaf, print, printAll, "V(E)EE", combi2E.size() * silhouetteEdges.size() * 2, 0, 2, 2, 0,
					e1, e2, true, combi2E, silhouetteEdges, silhouetteVertices)) return true;

	if (checkCombi(rst, prim, ray, leaf, print, printAll, "SEET", combi1S.size() * combi2E.size() * 3, 1, 0, 2, 1,
					splitLines, combi1S, true, combi2E, silhouetteEdges)) return true;

	if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSEE", combi2S.size() * combi2E.size(), 2, 0, 2, 0,
					splitLines, combi2S, true, combi2E, silhouetteEdges)) return true;
}

if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSET", combi2S.size() * combi1E.size() * 3, 2, 0, 1, 1,
				splitLines, combi2S, true, combi1E, silhouetteEdges)) return true;

std::vector<std::vector<int>> combi3S = Combinations::combi3(splitLines.size());
if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSSE", combi3S.size() * combi1E.size(), 3, 0, 1, 0,
				splitLines, combi3S, true, combi1E, silhouetteEdges)) return true;

return false;
}

bool RSTBuilderExact::check1Prim(RaySpaceTree* rst, Primitive* prim, Line4& ray, Node* leaf, bool print, int edgeSelection,
									bool getedges, std::vector<glm::vec3>& edges, std::vector<Ray>& esledges) {
	std::vector<Ray> splitLines;
	std::vector<bool> sideLines;

	rst->getSplittingLinesInLeafWithSide(leaf, splitLines, sideLines);
	if (!rst->alldir) {
		std::vector<Ray> boxSides = rst->model->boundingCube.getCubeSideSquare(rst->maindir).quadLines;
		for (int i = 0; i < 4; i++) {
			Ray boxside = boxSides[i];
			boxside.index = splitLines.size() + i;
			splitLines.push_back(boxside);
		}
	}

	std::vector<Ray> filteredsplitLines;
	std::vector<bool> filterdsideLines;
	filterSplittingLines(rst, leaf, splitLines, sideLines, filteredsplitLines, filterdsideLines);
	int size = filteredsplitLines.size();// +boxSides.size();
	if (size == 0) return false;
	std::vector<std::vector<int>> combi2 = Combinations::combi2(size);
	std::vector<std::vector<int>> combi3 = Combinations::combi3(size);
	std::vector<std::vector<int>> combi4 = Combinations::combi4(size);

	if (checkPrim(rst, prim, combi2, combi3, combi4, filteredsplitLines, filterdsideLines, ray, leaf, print, edgeSelection, getedges, edges, esledges)) return true;
	else if (print) std::cout << "Not Found" << std::endl;
	return false;
}

// Edgeselection 0 --> check from samples
// Edgeselection 1 --> check from parent node
bool RSTBuilderExact::checkLeaf(RaySpaceTree *rst, Node* node, std::vector<Line4>& rays, bool getrays, int edgeSelection, std::vector<int>& notfoundprim, bool print) {
	std::vector<Ray> splitLines;
	std::vector<bool> sideLines;

	rst->getSplittingLinesInLeafWithSide(node, splitLines, sideLines);
	if (!rst->alldir) {
		std::vector<Ray> boxSides = rst->model->boundingCube.getCubeSideSquare(rst->maindir).quadLines;
		for (int i = 0; i < 4; i++) {
			Ray boxside = boxSides[i];
			boxside.index = splitLines.size() + i;
			splitLines.push_back(boxside);
		}
	}
	std::vector<Ray> filteredsplitLines;
	std::vector<bool> filterdsideLines;
	filterSplittingLines(rst, node, splitLines, sideLines, filteredsplitLines, filterdsideLines);

	int size = filteredsplitLines.size();// +boxSides.size();
	std::vector<std::vector<int>> combi2 = Combinations::combi2(size);
	std::vector<std::vector<int>> combi3 = Combinations::combi3(size);
	std::vector<std::vector<int>> combi4 = Combinations::combi4(size);

	for (int i : node->primitiveSet) {
		Line4 ray;
		if (checkPrim(rst, rst->model->triangles[i], combi2, combi3, combi4, filteredsplitLines, filterdsideLines, ray, node, print, edgeSelection) && getrays) {
			rays.push_back(ray);
		}
		else {
			notfoundprim.push_back(i);
			std::cout << " DID NOT FIND PRIM " << i << " IN LEAF NR " << node->index << std::endl;
		}
	}
	if (notfoundprim.size() > 0) std::cout << " DID NOT FIND " << notfoundprim.size() << " OF " << node->primitiveSet.size() << " PRIMS IN LEAF NR " << node->index << std::endl;
	return notfoundprim.size() == 0;
}

void RSTBuilderExact::checkLeaves(RaySpaceTree* rst) {
	int leafcount = 0;
	for (Node* n : rst->nodes) {
		if (!n->leaf) continue;
		std::cout << "leaf index: " << n->index - rst->noLeaves << " from total: " << rst->noLeaves << std::endl;
		std::vector<Line4> rays;
		checkLeaf(rst, n, rays, false, 0);
	}
}

std::vector<Line4> RSTBuilderExact::getExtremalStabbingInLeaf(RaySpaceTree* rst, Node* n, std::vector<int>& notfoundprim, bool print) {
	std::vector<Line4> rays;
	checkLeaf(rst, n, rays, true, 0, notfoundprim, print);
	for (auto& r : rays) r.get3DfromPlucker();
	return rays;
}

bool RSTBuilderExact::checkEdgeSplittingDuplicates(Primitive* prim, std::vector<Ray>& splitLines, std::vector<bool>& sideLines) {
	for (int i = 0; i < sideLines.size(); i++) {
		for (Edge* e : prim->edges) {
			if (splitLines[i].equal(e->ray, 1E-5)) { // not very precise?!

				Ray checkRay = Ray(prim->center + prim->normal, prim->center);

				if (splitLines[i].side(checkRay) != sideLines[i]) {
					return false;
				}
			}
		}
	}
	return true;
}

bool RSTBuilderExact::checkPrim(RaySpaceTree* rst, Primitive* prim, std::vector<std::vector<int>>& combi2, std::vector<std::vector<int>>& combi3, std::vector<std::vector<int>>& combi4,
								std::vector<Ray> splitLines, std::vector<bool>& sideLines, Line4& ray, Node* leaf, bool print, int edgeSelection, bool getedges,
								std::vector<glm::vec3>& edges, std::vector<Ray>& eslEdges) {

	bool printAll = false;
	int splitLinesSize = splitLines.size();

	// Check if triangle edges are same as splitting lines and if yes, if it lies on the correct side
	if (!checkEdgeSplittingDuplicates(prim, splitLines, sideLines)) return false;

	// Check basic combis of extremal stabbing lines
	std::vector<Ray> e1;
	std::vector<std::vector<int>> e2;
	std::vector<Edge> e3;
	std::vector<int> e4;

	//if (!alldir && glm::dot(model->normalPerTri[prim], maindir) > 0) return false;

	// Check without occlusion to see if prim can be excluded
	if (!checkCombi(rst, prim, ray, leaf, false, printAll, "SSV(T)", combi2.size() * 3, 2, 0, 0, 2, splitLines, combi2, false) &&
		!checkCombi(rst, prim, ray, leaf, false, printAll, "SSST", combi3.size() * 3, 3, 0, 0, 1, splitLines, combi3, false) &&
		!checkCombi(rst, prim, ray, leaf, false, printAll, "SSSS", combi4.size(), 4, 0, 0, 0, splitLines, combi4, false)) return false;


	if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSV(T)", combi2.size() * 3, 2, 0, 0, 2, splitLines, combi2)) return true;
	if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSST", combi3.size() * 3, 3, 0, 0, 1, splitLines, combi3)) return true;

	// Find silhouette edges for primitive
	std::vector<Edge*> silhouetteEdges;
	std::vector<Edge*> silhouetteEdgesFirst;
	std::vector<Edge*> silhouetteEdgesSecond;

	if (edgeSelection == 0) {
		rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, silhouetteEdges);
		//exact would be to test if triangles belong in leaf
		for (Edge* e : silhouetteEdges) {
			bool found = false;
			for (Primitive* pr : e->triangles) {
				if (leaf->primitiveSet.find(pr->id) != leaf->primitiveSet.end()) {
					silhouetteEdgesFirst.push_back(e);
					found = true; break;
				}
			}
			if (!found) {
				Line4 eslEdge;
				if (checkEdgeInLeafCombis(rst, e, leaf, splitLines,combi2, eslEdge) || //bug????
					checkEdgeInLeafCombis(rst, e, leaf, splitLines, combi3, eslEdge)) {
					silhouetteEdgesSecond.push_back(e);
					if (getedges) {
						eslEdge.get3DfromPlucker();
						eslEdges.push_back(eslEdge);
					}
				}
			}
		}
	}
	else if (edgeSelection == 1)
		rst->model->findSilhouetteEdgesForTri(prim, rst->alldir, rst->maindir, silhouetteEdgesFirst, leaf->parent->primitiveSet);

	if (getedges) {
		for (Edge* e : silhouetteEdgesFirst) {
			for (Vertex* v : e->vertices) edges.push_back(v->pos);
		}
		for (Edge* e : silhouetteEdgesSecond) {
			for (Vertex* v : e->vertices) edges.push_back(v->pos);
		}
	}

	// Check all combis involving silhouette edges of some sort
	std::set<Vertex*, Vertex::cmp_ptr> silhEdgeVertices;
	silhouetteEdges = std::vector<Edge*>();

	if (silhouetteEdgesFirst.size() > 0) {
		if (checkSilhouetteCombis(rst, prim, ray, leaf, print, printAll, splitLines, silhouetteEdgesFirst, silhouetteEdges,
									silhEdgeVertices)) return true;
	}

	// Check basic (but large) combi of extremal stabbing lines
	if (checkCombi(rst, prim, ray, leaf, print, printAll, "SSSS", combi4.size(), 4, 0, 0, 0, splitLines, combi4)) return true;

	// Check combis involving second tier silhouette edges
	if (silhouetteEdgesSecond.size() > 0) {
		if (checkSilhouetteCombis(rst, prim, ray, leaf, print, printAll, splitLines, silhouetteEdgesSecond, silhouetteEdges,
									silhEdgeVertices)) return true;
	}

	return false;
}

