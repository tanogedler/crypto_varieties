#pragma once

#include <stdint.h>
#include <cstring>
#include <initializer_list>


union block512_t {
	uint8_t x8[64];
	uint64_t x64[8];
	block512_t(){}
	block512_t(const std::initializer_list<uint64_t> &x) {
		if (x.size() < 8) {
			memset(this, 0, sizeof(block512_t));
		}
		auto it = x.begin();
		size_t i = 0, size = x.size();
		while (i < size) {
			x64[i] = *it;
			++it; ++i;
		}
	}
	block512_t(uint8_t x) {
		for (int i = 0; i < 64; ++i) x8[i] = x;
	}
	bool operator==(const block512_t &x) const {
		for (int i = 0; i < 8; ++i) {
			if (x64[i] != x.x64[i]) {
				return false;
			}
		}
		return true;
	}
};

enum hash_mode_t
{
	hm512, hm256
};
