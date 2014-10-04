#pragma once

#include <cstdio>
#include "types.h"
#include "resource.h"

namespace hash_nm {
	block512_t hash(char *buf, size_t size, hash_mode_t mode = hm512);
	void precalc_mul_table();

//	#include "hash.hpp"
}
