
// CAS.cpp

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <ctime>
#include "polynomial.h"

using namespace std;
using namespace pol;


int main(int argc, char **argv) {
	if (argc > 0 && argv) {
		try {
			hash_nm::precalc_mul_table();
			FILE *fout = fopen("out", "wb");
			File in("in", READ);
			pol_t kx = get_rand_pol(KEY_DEG, false), ky = get_rand_pol(KEY_DEG, false);
			pol3v_t D = get_open_key(kx, ky);
			encrypt(in, fout, D);
			fclose(fout);
			in.~File();
			FILE *fin = fopen("out", "rb");
			File out("out2", WRITE);
			decrypt(fin, out, kx, ky);
			fclose(fin);
		}
		catch (char *str) {
			cout << str << endl;
			scanf("\n");
		}
	}
	return 0;
}
