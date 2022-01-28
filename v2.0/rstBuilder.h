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


	static Ray splitRandomEdge(RaySpaceTree* rst) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int tri = int(r * (rst->model->triangles.size()));
		float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		int edge = int(r2 * 3.f);
		return rst->model->triangles[tri]->edges[edge]->ray;
	}

	static Ray splitRandomVertexOrtho(RaySpaceTree* rst) {
		float r = (float)rand() / static_cast <float> (RAND_MAX);
		int r1 = int(r * rst->model->vertices.size());
		glm::vec3 pos = rst->model->vertices[r1]->pos;
		if (rand() % 2) return Ray(glm::vec3(pos.x, 0, pos.z), glm::vec3(pos.x, 1, pos.z));
		else return Ray(glm::vec3(pos.x, pos.y, 0), glm::vec3(pos.x, pos.y, 1));
	}

	//NEEDS FIXING
	static Ray splitRandomOrtho(RaySpaceTree* rst) {

		//float rX = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		//float rZ = -1.f + 2.f * float(rand()) / static_cast <float> (RAND_MAX);
		//float rY = 2.f * float(rand()) / static_cast <float> (RAND_MAX);

		//if (rand() % 2) return Ray(glm::vec3(rX, 0, rZ), glm::vec3(rX, 1, rZ));
		//else return Ray(glm::vec3(rX, rY, 0), glm::vec3(rX, rY, 1));
		return Ray();
	}

	static Ray splitRandomEdgeOffset(RaySpaceTree* rst) {
		//float r1 = (float)rand() / static_cast <float> (RAND_MAX);
		//int rtri = int(r1 * (rst->model->vertices.size() / 3));
		//float r2 = (float)rand() / static_cast <float> (RAND_MAX);
		//int redge = int(r2 * 3.f);
		//float r3 = (float)rand() / static_cast <float> (RAND_MAX);
		//float roffset = 0.1 * r3;
		//return Ray(rst->model->vertices[3 * rtri + redge].pos + roffset, rst->model->vertices[3 * rtri + (redge + 1) % 3].pos);
		return Ray();
	}

	static Ray getEdgeFromLeafNum(RaySpaceTree* rst, int leafNum) {
		return rst->model->triangles[leafNum / 3]->edges[leafNum % 3]->ray;
	}

	static void construct(RaySpaceTree* rst, int lvl, Node* node, Options::constructOption option, int splitnum) {
		if (lvl >= rst->depth) {
			node->leaf = true;
			return;
		}
		Ray splitter;
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
			case Options::FIRST_EDGES:
				splitter = getEdgeFromLeafNum(rst, node->index - rst->noLeaves);
		}

		bool dupli = false;
		// only bad if it's in the same subtree
		//for (Ray r : rst->splitters) if (splitter.equal(r, 1E-6)) { dupli = true; break; }
		if (dupli) construct(rst, lvl, node, option, splitnum);
		else {
			splitter.index = splitnum;
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

	static RaySpaceTree build(Model* model, int depth, bool alldir, char dir, int sgn, Options::BuildOptions& options, VisComponents& visComp) {
		glm::vec3 maindir = glm::vec3(0);
		if (!alldir) maindir = sgn * glm::ivec3(dir == 'X', dir == 'Y', dir == 'Z');

		model->findPotentialSilhouettes(alldir, maindir);

		RaySpaceTree rst(model, depth, alldir, maindir);
		if (options.construct != Options::ADAPTIVE) construct(&rst, options.construct);
		T::build(options, &rst, visComp);
		return rst;
	}

};
#endif