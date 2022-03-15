#ifndef RSTBUILDER_H
#define RSTBUILDER_H
#include "model.h"
#include "raytracerEmbree.h"
#include "buildOptions.h"
#include "visComponents.h"
#include "rst.h"

template <typename T>
class RSTBuilder {

public:
	//EmbreeTracer* embree;

	static void construct(RaySpaceTree* rst, Options::constructOption constructOption) {
		//int t = (int)time(NULL);
		//srand(t);
		std::vector<Ray> splitters;
		construct(rst, 0, rst->rootNode, constructOption, 1);
		rst->noLeaves = pow((int)2, rst->depth);
	};


	static Split splitRandomEdge(RaySpaceTree* rst) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int tri = int(r * (rst->model->triangles.size()));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int edge = int(r2 * 3.f);
		Edge* e = rst->model->triangles[tri]->edges[edge];
		return { e->ray, e };
	}

	static Split splitRandomVertexOrtho(RaySpaceTree* rst) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * rst->model->vertices.size());
		glm::vec3 pos = rst->model->vertices[r1]->pos;
		if (rand() % 2) return { Ray(glm::vec3(pos.x, 0, pos.z), glm::vec3(pos.x, 1, pos.z)) };
		else return { Ray(glm::vec3(pos.x, pos.y, 0), glm::vec3(pos.x, pos.y, 1)) };
	}

	//NEEDS FIXING
	static Split splitRandomOrtho(RaySpaceTree* rst) {

		//float rX = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		//float rZ = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		//float rY = 2.f * float(rand()) / static_cast <float> (RAND_MAX);

		//if (rand() % 2) return Ray(glm::vec3(rX, 0, rZ), glm::vec3(rX, 1, rZ));
		//else return Ray(glm::vec3(rX, rY, 0), glm::vec3(rX, rY, 1));
		return Split();
	}

	static Split splitRandomEdgeOffset(RaySpaceTree* rst) {
		//float r1 = (float)rand() / static_cast <float> (RAND_MAX);
		//int rtri = int(r1 * (rst->model->vertices.size() / 3));
		//float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		//int redge = int(r2 * 3.f);
		//float r3 = (float)rand() / static_cast <float> (RAND_MAX);
		//float roffset = 0.1 * r3;
		//return Ray(rst->model->vertices[3 * rtri + redge].pos + roffset, rst->model->vertices[3 * rtri + (redge + 1) % 3].pos);
		return Split();
	}

	static Split getEdgeFromLeafNum(RaySpaceTree* rst, int nodeNum) {
		Edge e* = rst->model->triangles[nodeNum / 3]->edges[nodeNum % 3];
		return { e->ray,  e};
	}

	static Split getEdgeFromLevel(RaySpaceTree* rst, int level) {
		int triStep = rst->model->edges.size() / rst->depth - 1;
		if (triStep == 0) std::cout << "not possible!" << std::endl;
		Edge* e = rst->model->edgeVector[level * triStep];
		return { e->ray, e };
	}

	static void construct(RaySpaceTree* rst, int lvl, Node* node, Options::constructOption option, int splitnum) {
		if (lvl >= rst->depth) {
			node->leaf = true;
			return;
		}
		Split splitter;
		switch(option) {
			case Options::RANDOM_EDGE:
				splitter = splitRandomEdge(rst);
				break;
			case Options::RANDOM_VERTEX_ORTHO:
				splitter = splitRandomVertexOrtho(rst);
				break;
			case Options::RANDOM_ORTHO:
				splitter = splitRandomOrtho(rst);
				break;
			case Options::SAME_LEVEL:
				splitter = getEdgeFromLevel(rst, lvl);
		}

		bool dupli = false;
		// only bad if it's in the same subtree
		std::vector<Split> parentSplitters;
		rst->getSplittingLinesInNode(node, parentSplitters);
		for (Split &ps : parentSplitters) if (splitter.ray.equal(ps.ray, 1E-6)) { dupli = true; break; }
		if (dupli) construct(rst, lvl, node, option, splitnum);
		else {
			splitter.id = splitnum;
			splitter.ray.index = splitnum;
			rst->splitters.push_back(splitter);

			node->splitter = splitter;
			Node* leftnode = new Node(node->index * 2 + 1, lvl + 1);
			node->leftNode = leftnode;
			leftnode->parent = node;
			rst->nodes.push_back(leftnode);
			construct(rst, lvl + 1, node->leftNode, option, splitnum * 2);

			Node* rightnode = new Node(node->index * 2 + 2, lvl + 1);
			rightnode->parent = node;
			node->rightNode = rightnode;
			rst->nodes.push_back(rightnode);
			construct(rst, lvl + 1, node->rightNode, option, splitnum * 2 + 1);
		}
	}

	static RaySpaceTree build(Model* model, char dir, int sgn, Options::BuildOptions& options, VisComponents& visComp) {
		glm::vec3 maindir = glm::vec3(0);
		if (!options.alldir) maindir = sgn * glm::ivec3(dir == 'X', dir == 'Y', dir == 'Z');

		model->findPotentialSilhouettes(options.alldir, maindir);

		RaySpaceTree rst(model, options.depth, options.alldir, maindir);
		if (options.construct != Options::ADAPTIVE) construct(&rst, options.construct);
		T::build(options, &rst, visComp);
		return rst;
	}

};
#endif