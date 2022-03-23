#pragma once

struct ExactBuildTracker {
	int totalTime = 0;
	int totalPrimitives = 0;
	int alreadyInNode = 0;
	int notInNodeFast = 0;
	int notInNodeSlow = 0;
	int inNodeFast = 0;
	int inNodeSlow = 0;
};