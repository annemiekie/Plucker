#ifndef COMBI_H
#define COMBI_H

#include <vector>

namespace Combinations {

	static std::vector<std::vector<int>> combi2(int size) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				std::vector<int> num = { i,j };
				combinations.push_back(num);
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi3(int size) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				for (int k = j + 1; k < size; k++) {
					std::vector<int> num = { i,j,k };
					combinations.push_back(num);
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi4(int size) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				for (int k = j + 1; k < size; k++) {
					for (int m = k + 1; m < size; m++) {
						std::vector<int> num = { i,j,k,m };
						combinations.push_back(num);
					}
				}
			}
		}
		return combinations;
	}

}
#endif