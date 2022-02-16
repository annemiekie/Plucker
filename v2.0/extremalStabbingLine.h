#pragma once
#include "lineThroughFourLines.h"
#include "node.h"
#include "primitive.h"

struct ExtremalStabbingLine {
	Node* leaf;
	Primitive* prim;
	Line4 ray;
};