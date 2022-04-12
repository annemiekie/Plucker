#include "rst.h"
#include <chrono>

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

void RaySpaceTree::putPrimitive(Ray& ray, int primId, bool putRay, bool putPrim, 
								bool putInNodes, int putFromDepth) {
	putPrimitive(ray, primId, rootNode, putRay, putPrim, putInNodes, putFromDepth);
}

void RaySpaceTree::putPrimitive(Ray& ray, int primId, Node* node, bool putRay, bool putPrim,
								bool putInNodes, int putFromDepth) {
	if (node->leaf) {
		if (putRay && putPrim) node->insert(primId, ray);
		else if (putPrim) node->insert(primId);
		else if (putRay) node->insert(ray);
		return;
	}
	else if (putInNodes && putPrim && node->depth >= putFromDepth) node->insert(primId);

	if (node->splitter.ray.side(ray)) 
		putPrimitive(ray, primId, node->leftNode, putRay, putPrim, putInNodes, putFromDepth);
	else putPrimitive(ray, primId, node->rightNode, putRay, putPrim, putInNodes, putFromDepth);
}

Node* RaySpaceTree::descend(Ray& ray)
{
	return descend(ray, rootNode);
}

Node* RaySpaceTree::descend(Ray& ray, Node* node)
{
	if (node->leaf) return node;
	if (node->splitter.ray.side(ray)) return descend(ray, node->leftNode);
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

void RaySpaceTree::getSplittingLinesInLeaf(int leafnum, std::vector<Split>& lines) {
	Node* leaf = getLeafFromNum(leafnum);
	getSplittingLinesInNode(leaf, lines);
};

void RaySpaceTree::getSplittingLinesInNode(Node* n, std::vector<Split>& lines) {
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
		lines.push_back(n->splitter.ray);
	}
	//getSplittingLines(rootNode, lines);
	return lines;
}

void RaySpaceTree::getSplittingLines(Node* node, std::vector<Ray>& lines) {
	if (!node->leaf) {
		lines.push_back(node->splitter.ray);
		getSplittingLines(node->leftNode, lines);
		getSplittingLines(node->rightNode, lines);
	}
}
//
//std::vector<glm::vec3> RaySpaceTree::getSplittingLinesInGeo(GeoObject* object) {
//
//	std::vector<Ray> lines = getAllSplittingLines();
//	std::vector<glm::vec3> splitters;
//	float start, end;
//
//	for (Ray line : lines) {
//		object->intersect(line, start, end);
//		splitters.push_back(line.origin + (double)start * line.direction);
//		splitters.push_back(line.origin + (double)end * line.direction);
//	}
//	return splitters;
//}

void RaySpaceTree::getSplittingLinesInNodeWithSide(Node* n, std::vector<SplitSide>& lines) {
	if (n == rootNode) return;
	else {
		SplitSide split(n->parent->splitter, (n->parent->leftNode == n));
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