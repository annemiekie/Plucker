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

	static std::vector<std::vector<int>> combi2(int size, int skip = 0) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				if (skip && j < skip) continue;
				std::vector<int> num = { i,j };
				combinations.push_back(num);
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi1(int size, int skip = 0) {
		std::vector<std::vector<int>> combinations;
		for (int i = skip; i < size; i++) {
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


	static std::vector<std::vector<int>> combiAddSelective(int size, std::vector<std::vector<int>>& combiChecked, int skip = 0) {
		std::vector<std::vector<int>> combinations;
		std::vector<std::vector<int>> combifind;
		int nrOfElInCombiM1 = combiChecked[0].size() - 1;
		int nrOfElInCombi = combiChecked[0].size();
		std::vector<int> check(nrOfElInCombi);
		
		if (combiChecked[0].size() == 2) combifind = combi1(2); 
		else if (combiChecked[0].size() == 3) combifind = combi2(3); 

		//for (std::vector<int>& cc : combiChecked) {
		for (int k = 0; k < combiChecked.size(); k++) {
			for (int i = combiChecked[k][nrOfElInCombiM1] + 1; i < size; i++) {
				if (i < skip) continue;
				int found = 0;

				for (std::vector<int>& cf : combifind) {

					for (int j = 0; j<combifind[0].size(); j++) check[j] = combiChecked[k][cf[j]];
					check[nrOfElInCombiM1] = i;

					for (int m = k + 1; m < combiChecked.size(); m++) {
						if (combiChecked[m] == check) {
							found++;
							break;
						}
					}
				}

				if (found == nrOfElInCombi) {
					std::vector<int> num;
					for (int c : combiChecked[k]) num.push_back(c);
					num.push_back(i);
					combinations.push_back(num);
				}
			}
		}
		return combinations;
	}

	static std::vector<std::vector<int>> combi3(int size, int skip = 1) {
		std::vector<std::vector<int>> combinations;
		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++) {
				for (int k = j + skip; k < size; k++) {
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