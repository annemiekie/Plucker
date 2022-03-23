#ifndef OPTIONS_H
#define OPTIONS_H
#include <unordered_map>

namespace Options {

	enum constructOption { ADAPTIVE, RANDOM_EDGE, RANDOM_ORTHO, RANDOM_VERTEX_ORTHO, FIRST_EDGES, SAME_LEVEL };
	static std::unordered_map<std::string, constructOption> const table = {
		{"ADAPTIVE",constructOption::ADAPTIVE},
		{"RANDOM_EDGE",constructOption::RANDOM_EDGE},
		{"RANDOM_ORTHO",constructOption::RANDOM_ORTHO},
		{"RANDOM_VERTEX_ORTHO",constructOption::RANDOM_VERTEX_ORTHO},
		{"FIRST_EDGES",constructOption::FIRST_EDGES},
		{"SAME_LEVEL",constructOption::SAME_LEVEL} };

	enum samplingType { UNIFORM_SPHERE, CUBE_DOMES };
	static std::unordered_map<std::string, samplingType> const table2 = {
		{"UNIFORM_SPHERE",samplingType::UNIFORM_SPHERE},
		{"CUBE_DOMES",samplingType::CUBE_DOMES} };

	struct BuildOptions {
		constructOption construct;
		samplingType samplingtype;
		int s_h;
		int s_w;
		int noSamples;
		int sampleStoreFillRatio;
		bool rasterizationSampling = false;
		bool storeSamples = false;
		bool storeAllSamples = false;
		bool cacheCombi = false;
		bool cacheEEE = false;
		bool cacheEE = false;
		bool fillMoreSamples = false;
		bool alldir = false;
		int depth = 0;
		bool sampling = false;
		bool exact = false;
		int exactStartLevel = 0;
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