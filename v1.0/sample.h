#pragma once
#include "ray.h"

struct Sample {
	Ray ray;
	int prim = -1;
};