#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <utility>
#include "resource.h"
#include "types.h"
#include "hash.h"

using namespace std;

uint64_t lpsbox[8][0x100];

void hash_nm::precalc_mul_table() {
	for (int i = 0; i < 8; ++i) {
		for (int j = 0; j < 0x100; ++j){
			uint64_t t = 0;
			uint8_t p = sbox[j];
			for (int k = 0; k < 8; ++k)
			if (p &(1 << k))
				t ^= lbox[0x3f - ((i << 3) | k)];
			lpsbox[i][j] = t;
		}
	}
}
void lps(block512_t &in) {
	block512_t res;
	for (int i = 0; i < 8; ++i){
		uint64_t t = lpsbox[0][in.x8[i]];
		for (int k = 1; k < 8; ++k)
			t ^= lpsbox[k][in.x8[i + k * 8]];
		res.x64[i] = t;
	}
	in = res;
}

void add512(block512_t &a, const block512_t &b) {
	uint32_t tmp = 0;
	for (int i = 0; i < 64; ++i) {
		tmp = a.x8[i] + b.x8[i] + (tmp >> 8);
		a.x8[i] = tmp & 0xFF;
	}
}
inline void x(block512_t &a, const block512_t &b) {
	a.x64[0] ^= b.x64[0]; a.x64[1] ^= b.x64[1]; a.x64[2] ^= b.x64[2]; a.x64[3] ^= b.x64[3];
	a.x64[4] ^= b.x64[4]; a.x64[5] ^= b.x64[5]; a.x64[6] ^= b.x64[6]; a.x64[7] ^= b.x64[7];
}
void g(const block512_t &n, block512_t &h, const block512_t &m) {
	block512_t res = h;
	x(res, n);
	lps(res);
	block512_t k = res;
	x(res, m);
	for (int i = 0; i < 12; ++i) {
		lps(res);
		x(k, c[i]);
		lps(k);
		x(res, k);
	}
	x(h, res);
	x(h, m);
}
block512_t hash_nm::hash(char *buf, size_t size, hash_mode_t mode) {
	block512_t 	m, h(mode == hm256), n(0), s(0);
	for (int64_t i = size - 64; i > -1; i -= 64) {
		memcpy(&m, buf + i, sizeof(m));
		g(n, h, m);
		add512(n, value512);
		add512(s, m);
	}
	memset(m.x8, 0, sizeof(m));
	memcpy(m.x8, buf + size, size % 64);
	m.x8[size % 64] = 1;
	g(n, h, m);
	add512(s, m);
	memset(m.x8, 0, sizeof(m));
	m.x64[0] = 8 * (size % 64);
	add512(n, m);
	m.x64[0] = 0;
	g(m, h, n);
	g(m, h, s);
	return h;
}
