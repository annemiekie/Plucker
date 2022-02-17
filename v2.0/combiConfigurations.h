#pragma once
#include "combinations.h"

struct CombiConfigurations {
	std::vector<std::vector<int>> c1, c2, c3, c4;

	CombiConfigurations() {};

	CombiConfigurations(int size) {
		c1 = Combinations::combi1(size);
		c2 = Combinations::combi2(size);
		c3 = Combinations::combi3(size);
		c4 = Combinations::combi4(size);
	};

};