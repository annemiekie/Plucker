#ifndef OPTIONS_H
#define OPTIONS_H
#include <unordered_map>

namespace Options {

	enum constructOption { ADAPTIVE, RANDOM_EDGE, RANDOM_ORTHO, RANDOM_VERTEX_ORTHO, RANDOM_EDGE_OFFSET };
	static std::unordered_map<std::string, constructOption> const table = {
		{"ADAPTIVE",constructOption::ADAPTIVE},
		{"RANDOM_EDGE",constructOption::RANDOM_EDGE},
		{"RANDOM_ORTHO",constructOption::RANDOM_ORTHO},
		{"RANDOM_VERTEX_ORTHO",constructOption::RANDOM_VERTEX_ORTHO},
		{"RANDOM_EDGE_OFFSET",constructOption::RANDOM_EDGE_OFFSET} };

	struct BuildOptions {
		constructOption construct;
		int s_h;
		int s_w;
		int noSamples;
		bool rasterizationSampling = false;
		bool storeRays = false;
		bool cacheCombi;
	};

	//constructOption findconstruct(std::string str) {
	//	auto it = table.find(str);
	//	if (it != table.end()) {
	//		return it->second;
	//	}
	//	else return constructOption::RANDOM_EDGE;
	//}
};

#endif