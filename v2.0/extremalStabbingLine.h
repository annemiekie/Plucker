#pragma once
#include "node.h"
#include "primitive.h"
#include "eslCandidate.h"

struct ExtremalStabbingLine : public ESLCandidate {
	std::vector<Primitive*> forPrimitives;
};