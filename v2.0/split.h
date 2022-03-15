#pragma once
#include "ray.h"
#include "edge.h"

struct Split {
	Ray ray;
	Edge* edge = NULL;
	int id;
};