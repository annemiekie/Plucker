#pragma once
#include "split.h"

struct SplitSide : Split {
	bool side;
	SplitSide(Split& s, bool side) : Split{ s.ray, s.edge, s.id }, side(side) {};
	SplitSide() {};
};