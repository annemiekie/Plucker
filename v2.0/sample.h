#pragma once
#include "ray.h"

struct Sample {
	Ray ray;
	int prim = -1;
};

struct SampleInd {
	int prim = -1;
	uint64_t raynr = -1;
};