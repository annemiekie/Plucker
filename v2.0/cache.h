#pragma once
#include <vector>
#include <set>
#include "RobinHood.h"
#include "ray.h"

template <class T>
class Cache {

public:
	robin_hood::unordered_map<uint64_t, T> cache;
	int hitcount = 0;
	std::set<size_t> hashes = std::set<size_t>();
	uint64_t lastKey = 0;

	uint64_t makeCombiKey(std::vector<uint64_t>& keyComponents) {
		uint64_t combiKey = 0;
		int bitshift = 64 / keyComponents.size();
		for (int i = 0; i < keyComponents.size(); i++) combiKey = combiKey | keyComponents[i] << i * bitshift;
		lastKey = combiKey;
		return combiKey;
	};

	bool getValue(std::vector<uint64_t>& keyComponents, T& value) {
		if (keyComponents.size() == 0) return false;
		uint64_t combiKey = makeCombiKey(keyComponents);
		if (!cache.contains(combiKey)) return false;
		value = cache[combiKey];
		hitcount++;
		return true;
	};

	void storeValueAtLastKey(T& value) {
		cache[lastKey] = value;
	}
};
