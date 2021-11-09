#ifndef COMBI_H
#define COMBI_H

#include <vector>

namespace Combinations {

	static std::vector<std::vector<int>> combi21(int size1, int size2) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size1; i++) {
			for (int j=i+1; j<size1; j++) {
				for (int k = size1; k < size1 + size2; k++) {
					std::vector<int> num = { i, j, k };
					combinations.push_back(num);
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi12(int size1, int size2) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size1; i++) {
			for (int j = size1; j < size1+size2; j++) {
				for (int k = j+1; k < size1 + size2; k++) {
					std::vector<int> num = { i, j, k };
					combinations.push_back(num);
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi13(int size1, int size2) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size1; i++) {
			for (int j = size1; j < size1 + size2; j++) {
				for (int k = j + 1; k < size1 + size2; k++) {
					for (int l = k + 1; l < size1 + size2; l++) {
						std::vector<int> num = { i, j, k, l };
						combinations.push_back(num);
					}
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi31(int size1, int size2) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size1; i++) {
			for (int j = i + 1; j < size1; j++) {
				for (int k = j + 1; k < size1; k++) {
					for (int l = size1; l < size1 + size2; l++) {
						std::vector<int> num = { i, j, k, l };
						combinations.push_back(num);
					}
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi22(int size1, int size2) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size1; i++) {
			for (int j = i + 1; j < size1; j++) {
				for (int k = size1; k < size1 + size2; k++) {
					for (int l = k+1; l < size1 + size2; l++) {
						std::vector<int> num = { i, j, k, l };
						combinations.push_back(num);
					}
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi11(int size1, int size2) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size1; i++) {
			for (int j = size1; j < size1+size2; j++) {
				std::vector<int> num = { i, j };
				combinations.push_back(num);
			}
		}
		return combinations;
	}

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

	static std::vector<std::vector<int>> combi1(int size) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size; i++) {
			std::vector<int> num = { i};
			combinations.push_back(num);
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combiAdd2(int size, std::vector<std::vector<int>>& combiChecked) {
		std::vector<std::vector<int>> combinations;

		for (std::vector<int>& cc : combiChecked) {
			for (int i = 0; i < size; i++) {
				for (int j = i + 1; j < size; j++) {
					std::vector<int> num;
					num.push_back(i);
					num.push_back(j);
					for (int c : cc) num.push_back(c+size);
					combinations.push_back(num);
				}
			}
		}
		return combinations;
	}


	static std::vector<std::vector<int>> combiAdd1(int size, std::vector<std::vector<int>>& combiChecked) {
		std::vector<std::vector<int>> combinations;

		for (std::vector<int>& cc : combiChecked) {
			for (int i = 0; i < size; i++) {
				std::vector<int> num;
				num.push_back(i);
				for (int c : cc) num.push_back(c+size);
				combinations.push_back(num);
			}
		}
		return combinations;
	}


	static std::vector<std::vector<int>> combiAddSelective(int size, std::vector<std::vector<int>>& combiChecked) {
		std::vector<std::vector<int>> combinations;
		std::vector<std::vector<int>> combifind;
		if (combiChecked[0].size() == 2) combifind = combi1(2); 
		else if (combiChecked[0].size() == 3) combifind = combi2(3); 

		for (std::vector<int>& cc : combiChecked) {
			for (int i = cc[cc.size() - 1] + 1; i < size; i++) {
				int found = 0;
				for (std::vector<int>& cf : combifind) {
					bool add = false;
					std::vector<int> check;
					for (int f : cf) check.push_back(cc[f]); 
					check.push_back(i);  
					for (std::vector<int>& cc2 : combiChecked) {
						if (cc2 == check) {
							found++;
							break;
						}
					}
				}

				if (found == cc.size()) {
					std::vector<int> num;
					for (int c : cc) num.push_back(c);
					num.push_back(i);
					combinations.push_back(num);
				}
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