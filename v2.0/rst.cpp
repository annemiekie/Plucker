#include "rst.h"
#include <chrono>


//void RaySpaceTree::fillExact() {
//	Ray r;
//	//if (alldir) {
//	//#pragma omp parallel
//	for (Node* node : nodes) { //int n = 0; n < nodes.size(); n++) {// 
//		std::cout << "Computing node nr: " << node->index << " of " << nodes.size() << std::endl;
//		if (node->leaf) {
//			//std::cout << "Computing node nr: " << nodes[n]->index << " of " << nodes.size() << std::endl;
//			for (int i = 0; i < model->primsize; i++) {
//				std::cout << i << ", ";
//				if (node->primitiveSet.find(i) == node->primitiveSet.end()) {
//					//if (glm::dot(model->normalPerTri[i], maindir) >= 0) continue;
//					//std::cout << "Primitive: " << node->index << " of " << model->primsize << ": " << std::endl;
//					if (check1Prim(i, r, node, false, 0))  node->primitiveSet.insert(i);
//				}
//			}
//
//		}
//		std::cout << std::endl;
//		if (cacheCombi)  std::cout << "Combi Cache, Size: " << combiCache.cache.size() << " Hits: " << combiCache.hitcount << " Hash Hit Size: " << combiCache.hashes.size() << std::endl;
//		if (model->cacheEEE) std::cout << "Edge Edge Edge Cache, Size: " << model->edgeEdgeEdgeCombis.size() << " Hits: " << model->cachehiteee << std::endl;
//		if (model->cacheEE) std::cout << "Edge Edge Cache, Size: " << model->edgeEdgeCombis[0].size() << " Hits: " << model->cachehitee << std::endl;
//	}
	//}
	//else {
	//	std::cout << "Computing node nr: 0 of " << nodes.size() << std::endl;
	//	std::cout << "Prim size: " << model->primsize << std::endl;
	//	//#pragma omp parallel
	//	for (int i = 0; i < model->primsize; i++) {
	//		//if (!alldir) 
	//		//if (glm::dot(model->normalPerTri[i], maindir)) continue;
	//		std::cout << i << ", ";
	//		if (check1Prim(i, r, rootNode, false, 0)) rootNode->primitiveSet.insert(i);//, 
	//		//else  rootNode->primitiveSet.insert(i);
	//	}
	//	std::cout << std::endl;
	//	fillExact(rootNode->leftNode);
	//	fillExact(rootNode->rightNode);
	//}

////
//void RaySpaceTree::fillExact(Node* node) {
//
//	std::cout << "Computing node nr: " << node->index << " of " << nodes.size() << std::endl;
//	std::cout << "Prim size: " << node->parent->primitiveSet.size() << std::endl;
//	Ray r;
//	//#pragma omp parallel
//	int cnt = 0;
//	for (int i : node->parent->primitiveSet) {
//		std::cout << cnt << ", ";
//		if (check1Prim(i, r, node, false, 1)) node->primitiveSet.insert(i);//, 1
//		cnt++;
//	}
//	std::cout << std::endl;
//
//	if (!node->leaf) {
//		fillExact(node->leftNode);
//		fillExact(node->rightNode);
//	}
//}




Ray RaySpaceTree::getRay(std::vector<Camera*>& cams, int raynr, glm::ivec2 res) {
	int camNr = raynr / (res.x * res.y);

	int rayInCam = raynr % (res.x * res.y);
	int y = rayInCam / res.x;
	int x = rayInCam % res.x;
	const glm::vec2 pixpos{ (x + .5f) / res.x * 2.0f - 1.0f, 1.0f - (y + .5f) / res.y * 2.0f };
	return cams[camNr]->pixRayDirection(pixpos);
}

bool RaySpaceTree::intoLeftNode(Ray& splitter, std::vector<Camera*>& cams, int raynr, glm::ivec2 res) {
	Ray ray = getRay(cams, raynr, res);
	return splitter.side(ray);
}


void RaySpaceTree::printTree() {
	printTree(rootNode, 0);
	std::cout << std::endl << std::endl;
}

void RaySpaceTree::printTree(Node* node, int level) {
	if (!node->leaf) {
		level++;
		std::cout << std::endl;
		for (int x = 0; x < level; x++) std::cout << "  ";
		std::cout << "L";
		printTree(node->leftNode, level);
		std::cout << std::endl;
		for (int x = 0; x < level; x++) std::cout << "  ";
		std::cout << "R";
		printTree(node->rightNode, level);
	}
	else {
		std::cout << " Leaf: tri = " << node->primitiveSet.size() << ", samples = " << node->primAndRayVector.size();
	}
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, bool putRay, bool putPrim) {
	putPrimitive(ray, primId, rootNode, putRay, putPrim);
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, Node* node, bool putRay, bool putPrim) {
	if (node->leaf) {
		if (putRay && putPrim) node->insert(primId, ray);
		else if (putPrim) node->insert(primId);
		else if (putRay) node->insert(ray);
		return;
	}
	if (node->splitter.side(ray)) putPrimitive(ray, primId, node->leftNode, putRay, putPrim);
	else putPrimitive(ray, primId, node->rightNode, putRay, putPrim);
}

Node* RaySpaceTree::descend(Ray& ray)
{
	return descend(ray, rootNode);
}

Node* RaySpaceTree::descend(Ray& ray, Node* node)
{
	if (node->leaf) return node;
	if (node->splitter.side(ray)) return descend(ray, node->leftNode);
	else return descend(ray, node->rightNode);
}

void RaySpaceTree::printLeafNodes() {
	for (Node* n : nodes) {
		if (n->leaf) {
			std::cout << n->index << ": ";
			for (int i : n->primitiveSet) {
				std::cout << i << ", ";
			}
			std::cout << std::endl;
		}
	}
}

std::vector<int> RaySpaceTree::countDuplicates(int size, std::vector<int>& nodenr) {
	std::vector<int> duplicates = std::vector<int>(size);
	countDuplicates(rootNode, duplicates, nodenr);
	return duplicates;
}

void RaySpaceTree::countDuplicates(Node* node, std::vector<int>& duplicates, std::vector<int>& nodenr) {
	if (!node->leaf) {
		countDuplicates(node->leftNode, duplicates, nodenr);
		countDuplicates(node->rightNode, duplicates, nodenr);
	}
	else {
		for (int x : node->primitiveSet) {
			duplicates[x] += 1;
			nodenr[node->index] += 1;
		}
	}
}

int RaySpaceTree::getNumberOfTriInleaf(int leafnum) {
	Node* n = getLeafFromNum(leafnum);
	return n->primitiveSet.size();
}

void RaySpaceTree::getSplittingLinesInLeaf(int leafnum, std::vector<Ray>& lines) {
	Node* leaf = getLeafFromNum(leafnum);
	getSplittingLinesInNode(leaf, lines);
};

void RaySpaceTree::getSplittingLinesInNode(Node* n, std::vector<Ray>& lines) {
	if (n == rootNode) return;
	else {
		lines.push_back(n->parent->splitter);
		return getSplittingLinesInNode(n->parent, lines);
	}
};

std::vector<Ray> RaySpaceTree::getAllSplittingLines() {
	std::vector<Ray> lines = std::vector<Ray>();
	for (Node* n : nodes) {
		if (n->leaf) continue;
		lines.push_back(n->splitter);
	}
	//getSplittingLines(rootNode, lines);
	return lines;
}

void RaySpaceTree::getSplittingLines(Node* node, std::vector<Ray>& lines) {
	if (!node->leaf) {
		lines.push_back(node->splitter);
		getSplittingLines(node->leftNode, lines);
		getSplittingLines(node->rightNode, lines);
	}
}

std::vector<glm::vec3> RaySpaceTree::getSplittingLinesInGeo(GeoObject* object) {

	std::vector<Ray> lines = getAllSplittingLines();
	std::vector<glm::vec3> splitters;
	float start, end;

	for (Ray line : lines) {
		object->intersect(line, start, end);
		splitters.push_back(line.origin + (double)start * line.direction);
		splitters.push_back(line.origin + (double)end * line.direction);
	}
	return splitters;
}

void RaySpaceTree::getSplittingLinesInNodeWithSide(Node* n, std::vector<Split>& lines) {
	if (n == rootNode) return;
	else {
		Split split = { n->parent->splitter, (n->parent->leftNode == n) };
		lines.push_back(split);
		return getSplittingLinesInNodeWithSide(n->parent, lines);
	}
}



int RaySpaceTree::numOfLeaf(int index) {
	int count = 0;
	for (Node* n : nodes) {
		if (n->leaf) {
			count++;
			if (n->index == index) return count;
		}
	}
	return 0;
}

Node* RaySpaceTree::getLeafFromNum(int leafNum) {
	int leaf = 0;
	for (Node* n : nodes) {
		if (n->leaf) {
			if (leaf == leafNum) {
				return n;
			}
			leaf++;
		}
	}
	return NULL;
}

std::vector<Ray> RaySpaceTree::getViewingLinesInLeaf(int leafNum)
{
	Node* node = getLeafFromNum(leafNum);
	return getViewingLinesInLeaf(node);
}

std::vector<Ray> RaySpaceTree::getViewingLinesInLeaf(Node* node)
{
	std::vector<Ray> rays;
	for (Sample& p : node->primAndRayVector) rays.push_back(p.ray);
	return rays;
}

std::vector<Ray> RaySpaceTree::getviewingLinesInLeafInTri(Node* node, int prim)
{
	std::vector<Ray> rays;
	for (Sample& p : node->primAndRayVector) {
		if (p.prim == prim) rays.push_back(p.ray);
	}
	return rays;
}