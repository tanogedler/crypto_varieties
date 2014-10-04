#include "polynomial.h"
//using namespace pol;
template <unsigned KEY_DEG>

void encrypt(File &in, FILE *out, const pol3v_t &D) {
	std::vector<u61_t> m;
	u61_t buf;
	u64_t size = 0, block_size;
	while (block_size = in.read(buf), block_size) {
		m.push_back(buf);
		size += block_size;
	}
	block512_t h = hash_nm::hash((char *)(m.data()), m.size() * sizeof(u61_t));
	fwrite(&h, sizeof(h), 1, out);
	u64_t max_deg = (m.size() - 1) / KEY_DEG + 2;
	pol3v_t f = XXX_get_f(max_deg, KEY_DEG);
	pol3v_t mes = XXX_get_mes(max_deg, KEY_DEG, size, m);
	(mes += f * D.get_same_rand_pol() + D * f.get_same_rand_pol()).write(out, TEXTPR);
	(mes += f * D.get_same_rand_pol() + D * f.get_same_rand_pol()).write(out, TEXTPR);
}

template <unsigned CLOSE_KEY_DEG, unsigned OPEN_KEY_DEG>
void decrypt(FILE *in, File &out, const pol_t &x, const pol_t &y) {
	block512_t h1, h2;
	fread(&h1, sizeof(h1), 1, in);
	pol3v_t tmp;
	tmp.read(in, TEXTPR);
	pol_t m, m1, fc, pol_x = get_x(), kxy = x * y;
	m1 = tmp.subst(x, y);
	tmp.read(in, TEXTPR);
	fc = tmp.subst(x, y) - m1;
	vector_t factor;
	fc.fact(CLOSE_KEY_DEG * OPEN_KEY_DEG * 2 + 1, factor);
	XXX_get_fact(factor, CLOSE_KEY_DEG * OPEN_KEY_DEG * 2 + 1, factor);
	std::vector<u61_t> res;
	u61_t a = x[CLOSE_KEY_DEG].inverse() * y[CLOSE_KEY_DEG].inverse();
	for (u32_t j = 0; j < factor.size(); ++j) {
		res.clear();
		m = m1 % (fc / factor[j]);
		u32_t deg = m.get_deg();
		u64_t size = m[deg] * a.pow(deg / (2 * CLOSE_KEY_DEG));
		m %= kxy.pow(deg / (2 * CLOSE_KEY_DEG)) * pol_x.pow(deg % (2 * CLOSE_KEY_DEG));
		deg = m.get_deg();
		for (u32_t i = 0; i <= deg; ++i) {
			res.push_back(m[deg - i] * a.pow((deg - i) / (2 * CLOSE_KEY_DEG)));
			m %= kxy.pow((deg - i) / (2 * CLOSE_KEY_DEG)) * pol_x.pow((deg - i) % (2 * CLOSE_KEY_DEG));
		}
		block512_t h2 = hash_nm::hash((char *)(res.data()), res.size() * sizeof(u61_t));
		if (h1 == h2) {
//			for (u32_t i = 0; i < res.size(); ++i) {
//				if (size > MODSIZE - 1) {
//					size -= out.write(res[i]);
//				}
//				else {
//					size -= out.write(res[i], size);
//				}
//			}
			out.flush();
			break;
		}
	}
}

template <unsigned CLOSE_KEY_DEG, unsigned OPEN_KEY_DEG>
   pol3v_t get_open_key(const pol_t &x, const pol_t &y) {
	polmap_t tmp;
	pol3v_t tmp2;
	for (u32_t i = 1; i <= OPEN_KEY_DEG; ++i) {
		tmp.insert(pol3v1_t(deg_t(i, i), get_rand_pol(2 * (OPEN_KEY_DEG - i) * CLOSE_KEY_DEG + 1, false)));
	}
	tmp.insert(pol3v1_t(deg_t(OPEN_KEY_DEG - 1, OPEN_KEY_DEG), get_rand_pol(1 + CLOSE_KEY_DEG, false)));
	tmp.insert(pol3v1_t(deg_t(OPEN_KEY_DEG, OPEN_KEY_DEG - 1), get_rand_pol(1 + CLOSE_KEY_DEG, false)));
	tmp.insert(pol3v1_t(deg_t(0, 1), get_rand_pol(CLOSE_KEY_DEG * (2 * OPEN_KEY_DEG - 1) + 1, false)));
	tmp.insert(pol3v1_t(deg_t(1, 0), get_rand_pol(CLOSE_KEY_DEG * (2 * OPEN_KEY_DEG - 1) + 1, false)));
	tmp2 = tmp;
	tmp.insert(pol3v1_t(deg_t(0, 0), pol_t() - tmp2.subst(x, y)));
	return tmp;
}
