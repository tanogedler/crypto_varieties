// polynomial.cpp

#define _CRT_SECURE_NO_WARNINGS

#include "polynomial.h"
#include <array>
#include <map>
#include <cstdio>

using namespace pol;

u61_t::u61_t(u64_t x):
	data(x % MOD) {
}
inline u61_t::operator pol::u64_t() const {
	return data;
}
inline u61_t u61_t::operator+(const u61_t &x) const {
	return (data + x.data) % MOD;
}
inline u61_t u61_t::operator-(const u61_t &x) const {
	return (data + MOD - x.data) % MOD;
}
u61_t u61_t::operator*(const u61_t &x) const {
//#pragma warning(push)
//#pragma warning(disable:4127)
	if (BUFSIZE > 4) {
//#pragma warning(pop)
		u64_t h1 = data >> 32, h2 = x.data >> 32, l1 = data & M32, l2 = x.data & M32;
		u64_t bot = l1 * l2;
		u64_t mid = h1 * l2 + h2 * l1 + (bot >> 32);
		bot = (mid << 32) + (bot & M32);
		u64_t top = ((h1 * h2 + (mid >> 32)) << (64 - MODSIZE)) + (bot >> MODSIZE);
		return (top + (bot & MOD)) % MOD;
	} else {
		return data * x.data % MOD;
	}
}
u61_t u61_t::operator/(const u61_t &x) const {
	return *this * (x.inverse());
}
u61_t u61_t::inverse() const {
	if (!data) {
		throw "Invalid inverse\n";
	}
	u64_t a = MOD, b = data, r;
	std::array<u64_t, 64> q;
	q[0] = 1;
	u32_t i = 0;
	do {
		++i;
		r = a % b;
		q[i] = a / b;
		a = b; b = r;
	} while (r);
	for (u32_t j = 2; j < i; ++j) {
		q[j] = q[j] * q[j - 1] + q[j - 2];
	}

	if (i & 1) {
		return u61_t(q[i - 1]);
	} else {
		return u61_t(MOD - q[i - 1]);
	}
}
u61_t& u61_t::operator+=(const u61_t &x) {
	return *this = *this + x;
}
u61_t& u61_t::operator-=(const u61_t &x) {
	return *this = *this - x;
}
u61_t& u61_t::operator*=(const u61_t &x) {
	return *this = *this * x;
}
u61_t& u61_t::operator/=(const u61_t &x) {
	return *this = *this / x;
}
inline u64_t* u61_t::get_addr() {
	return &data;
}
inline const u64_t* u61_t::get_addr() const {
	return &data;
}
inline u64_t& u61_t::get_buf() {
	return data;
}
inline const u64_t& u61_t::get_buf() const {
	return data;
}
void u61_t::print(FILE *fout, u32_t i) const {
	if ((i > 0 && data != 1) || !i) {
		fprintf(fout, "%llu", data);
	}
	if (i > 1) {
		fprintf(fout, "t<sup>%u</sup>", i);
	} else if (i == 1) {
		fprintf(fout, "t");
	}
}
u61_t u61_t::pow(u32_t p) const {
	u61_t res, tmp = *this;
	res = 1;
	while (p > 0) {
		if (p & 1) {
			res *= tmp;
		}
		p /= 2;
		tmp *= tmp;
	}
	return res;
}

pol_t::pol_t():
	data(1),
	degree(0)
{
	data[0] = 0;
}
void pol_t::change_degree(u32_t ndegree) {
	degree = ndegree;
	data.resize(degree + 1);
	for (u32_t i = 0; i <= degree; ++i) {
		data[i] = 0;
	}
}
pol_t& pol_t::operator=(const pol_t &x) {
	if (&x != this) {
		change_degree(x.degree);
		for (u32_t i = 0; i <= degree; ++i) {
			data[i] = x[i];
		}
	}
	return *this;
}
u61_t& pol_t::operator[](u32_t i) {
	if (i > degree) {
		throw "Out of array bound";
	}
	return data[i];
}
const u61_t& pol_t::operator[](u32_t i) const {
	if (i > degree) {
		throw "Out of array bound";
	}
	return data[i];
}
u64_t& pol_t::operator()(u32_t i) {
	if (i > degree) {
		throw "Out of array bound";
	}
	return data[i].get_buf();
}
const u64_t& pol_t::operator()(u32_t i) const {
	if (i > degree) {
		throw "Out of array bound";
	}
	return data[i].get_buf();
}
inline u32_t pol_t::get_deg() const {
	return degree;
}
pol_t pol_t::operator+(const pol_t &x) const{
	pol_t tmp;
	u32_t tempdeg;
	if (degree < x.degree) {
		tempdeg = x.degree;
	} else {
		tempdeg = degree;
	}
	tmp.change_degree(tempdeg);
	for (u32_t i = 0; i <= degree; ++i) {
		tmp[i] = data[i];
	}
	for (u32_t i = 0; i <= x.degree; ++i) {
		tmp[i] += x[i];
	}
	while (tmp.degree > 0 && !tmp[tmp.degree]) {
		--(tmp.degree);
	}
	return tmp;
}
pol_t pol_t::operator-(const pol_t &x) const {
	pol_t tmp;
	u32_t tempdeg;
	if (degree < x.degree) {
		tempdeg = x.degree;
	} else {
		tempdeg = degree;
	}
	tmp.change_degree(tempdeg);
	for (u32_t i = 0; i <= degree; ++i) {
		tmp[i] = data[i];
	}
	for (u32_t i = 0; i <= x.degree; ++i) {
		tmp[i] -= x[i];
	}
	while (tmp.degree > 0 && !tmp[tmp.degree]) {
		--(tmp.degree);
	}
	return tmp;
}
pol_t pol_t::operator*(const pol_t &x) const {
	pol_t tmp;
	tmp.change_degree(degree + x.degree);
	for (u32_t i = 0; i <= degree; ++i) {
		for (u32_t j = 0; j <= x.degree; ++j) {
			tmp[i + j] += data[i] * x[j];
		}
	}
	while (tmp.degree > 0 && !tmp[tmp.degree]) {
		--tmp.degree;
	}
	return tmp;
}
pol_t pol_t::operator%(const pol_t &x) const {
	if (!x) {
		throw "Polinomial denominator = 0";
	}
	pol_t tmp = *this;
	if (!x.degree) {
		return pol_t();
	} else {
		u61_t inv = x[x.degree].inverse();
		for (u32_t i = tmp.degree; i >= x.degree; --i) {
			u61_t q = tmp[i] * inv;
			for (u32_t j = 0; j <= x.degree; ++j) {
				tmp[i - j] -= x[x.degree - j] * q;
			}
		}
		while (tmp.degree > 0 && !tmp[tmp.degree]) {
			--tmp.degree;
		}
		return tmp;
	}
}
pol_t pol_t::operator/(const pol_t &x) const {
	if (!x) {
		throw "Polinomial devision by 0";
	} else if (x.degree > degree) {
		return pol_t();
	} else if (x.degree) {
		pol_t tmp = *this, res;
		res.change_degree(tmp.degree - x.degree);
		u61_t q = x[x.degree].inverse();
		for (u32_t i = tmp.degree; i >= x.degree; --i) {
			res[i - x.degree] = tmp[i] * q;
			for (u32_t j = 0; j <= x.degree; ++j) {
				tmp[i - j] -= x[x.degree - j] * res[i - x.degree];
			}
		}
		return res;
	} else {
		pol_t res = *this;
		u61_t q = x[0].inverse();
		for (u32_t i = 0; i <= res.degree; ++i) {
			res[i] *= q;
		}
		return res;
	}
}
pol_t& pol_t::operator+=(const pol_t&x) {
	return *this = *this + x;
}
pol_t& pol_t::operator-=(const pol_t&x) {
	return *this = *this - x;
}
pol_t& pol_t::operator*=(const pol_t&x) {
	return *this = *this * x;
}
pol_t& pol_t::operator%=(const pol_t&x) {
	return *this = *this % x;
}
pol_t& pol_t::operator/=(const pol_t&x) {
	return *this = *this / x;
}
pol_t::operator bool() const {
	return degree || data[0];
}
bool pol_t::operator==(const pol_t &x) const {
	if (x.degree != degree) {
		return false;
	} else {
		for (u32_t i = 0; i <= degree; ++i) {
			if (data[i] != x(i)) {
				return false;
			}
		}
	}
	return true;
}
pol_t pol_t::pow(u32_t p) const {
	pol_t res, tmp = *this;
	res(0) = 1;
	while (p > 0) {
		if (p & 1) {
			res *= tmp;
		}
		p /= 2;
		tmp *= tmp;
	}
	return res;
}
pol_t pol_t::mpow(u64_t p, const pol_t &m) const {
	pol_t res, tmp = *this % m;
	res(0) = 1;
	while (p > 0) {
		if (p & 1) {
			res = res * tmp % m;
		}
		p /= 2;
		tmp = tmp * tmp % m;
	}
	return res;
}
pol_t pol_t::mpow(u64_t p, u64_t d, const pol_t &m) const {
	pol_t tmp = mpow((p - 1) / 2, m), res;
	res[0] = 1;
	for (u64_t i = 0; i < d; ++i) {
		res = res * tmp % m;
		tmp = tmp.mpow(p, m);
	}
	return res;
}
pol_t pol_t::normalize() const {
	pol_t res = *this;
	u61_t q = data[degree].inverse();
	for (u32_t i = 0; i <= degree; ++i) {
		res[i] *= q;
	}
	return res;
}
pol_t pol_t::deriv() const{
	pol_t res;
	if (degree) {
		res.change_degree(degree - 1);
		for (u32_t i = degree; i > 0; --i) {
			res[i - 1] = u61_t(i) * data[i];
		}
	}
	while (res.degree > 0 && !res[res.degree]) {
		--res.degree;
	}
	return res;
}

void pol_t::read(FILE *fin, pol_format_t format) {
	degree = 0;
	switch (format)
	{
	case BINPR:
		std::fread(&degree, 1, 1, fin);
		if (degree > 256) {
			throw "Incorrect polynomial degree";
		}
		data.resize(degree + 1);
		for (u32_t i = 0; i <= degree; ++i) {
			u64_t buf;
			std::fread(&buf, BUFSIZE, 1, fin);
			if (buf >= MOD) {
				throw "Incorrect polinomial coefficient";
			}
			data[i] = buf;
		}
		break;
	case TEXTPR:
		std::fscanf(fin, "%u", &degree);
		if (degree > 256) {
			throw "Incorrect polynomial degree";
		}
		data.resize(degree + 1);
		for (u32_t i = 0; i <= degree; ++i) {
			u64_t buf;
			std::fscanf(fin, "%llu", &buf);
			if (buf >= MOD) {
				throw "Incorrect polinomial coefficient";
			}
			data[i] = buf;
		}
		break;
	case HTMLPR:
		break;
	default:
		break;
	}
}
void pol_t::write(FILE *fout, pol_format_t format, bool eol) const {
	switch (format)
	{
	case pol::BINPR:
		fwrite(&degree, 1, 1, fout);
		for (u32_t i = 0; i <= degree; ++i) {
			fwrite(data[i].get_addr(), BUFSIZE, 1, fout);
		}
		break;
	case pol::TEXTPR:
		fprintf(fout, "%u ", degree);
		for (u32_t i = 0; i <= degree; ++i) {
			fprintf(fout, "%llu ", data[i].get_buf());
		}
		if (eol) fprintf(fout, "\n");
		break;
	case pol::HTMLPR:
		{
			u32_t i = degree;
			u32_t min = 0;
			while (min < degree && !data[min]) ++min;
			--min;
			data[i].print(fout, i);
			--i;
			for (; i != min; --i) {
				if (data[i]) {
					fprintf(fout, " + ");
					data[i].print(fout, i);
				}
			}
			if (eol) fprintf(fout, "<br>");
		}
		break;
	default:
		break;
	}
}
void pol_t::ddf(u32_t max_deg, vector_t &res) const {
	pol_t v = *this, x = get_x(), xp = x, tmp;
	res.resize(max_deg + 1);
	for (u32_t i = 0; i <= max_deg; ++i) {
		res[i](0) = 1;
	}
	u32_t i = 1;
	for (; i <= max_deg; ++i)  {
		xp = xp.mpow(MOD, v);
		tmp = gcd(xp - x, v);
		if (tmp.degree > 0) {
			res[i] = tmp;
			v /= tmp;
		}
		fprintf(stderr, "%d of %d\n", i, max_deg);
	}
	res[0] = v;
}
void pol_t::split(u32_t deg, vector_t &res) const {
	u32_t count = degree / deg, found = 0;
	pol_t r_pol, r_gcd, unit;
	unit[0] = 1;
	std::list<pol_t> l1, l2;
	l1.push_back(*this);
	std::list<pol_t> *cur = &l1, *next = &l2;
	while (found < count) {
		r_pol = get_rand_pol(deg);
		next->clear();
		for (auto i = cur->begin(), end = cur->end(); i != end; ++i) {
			if (i->degree == deg) {
				res.push_back(*i);
				++found;
			} else {
				r_gcd = gcd(r_pol.mpow(MOD, deg, *i) - unit, *i);
				if (r_gcd.degree && r_gcd.degree < i->degree) {
					*i /= r_gcd;
					if (i->degree == deg) {
						res.push_back(*i);
						++found;
					}
					else {
						next->push_back(*i);
					}
					if (r_gcd.degree == deg) {
						res.push_back(r_gcd);
						++found;
					}
					else {
						next->push_back(r_gcd);
					}
				}
				else {
					next->push_back(*i);
				}
			}
		}
		std::swap(cur, next);
	}
}
void pol_t::fact(u32_t max_deg, vector_t &res) const {
	pol_t d = gcd(deriv(), *this), tmp;
	if (d.degree) {
		tmp = *this / d;
	} else {
		tmp = *this;
	}
	std::vector<pol_t> ddf_pol;
	tmp.ddf(max_deg, ddf_pol);
	for (u32_t i = 1; i <= max_deg; ++i) {
		if (ddf_pol[i].degree > 0) {
			ddf_pol[i].split(i, res);
		}
	}
	if (d.degree) {
		for (size_t i = 0, size = res.size(); i < size && d.degree; ++i) {
			while (!(d % res[i])) {
				res.push_back(res[i]);
				d /= res[i];
			}
		}
	}
}

deg_t::deg_t() : key(0) {}
deg_t::deg_t(u32_t a, u32_t b) :
	x(a),
	y(b)
{
}
bool deg_t::operator<(const deg_t &deg) const{
	return key < deg.key;
}
bool deg_t::operator==(const deg_t &deg) const{
	return key == deg.key;
}
deg_t deg_t::operator+(const deg_t &deg) const {
	return deg_t(x + deg.x, y + deg.y);
}

pol3v_t::pol3v_t() {
	data.clear();
	data.push_back(pol3v1_t());
	dx = 0; dy = 0; dt = 0;
}
pol3v_t::pol3v_t(const polmap_t &m) :
	data(),
	dx(0),
	dy(0),
	dt(0)
{
	data.clear();
	bool empty = true;
	for (auto i = m.begin(); i != m.end(); ++i) {
		if (i->second) {
			if (empty) {
				empty = false;
			}
			data.push_back(*i);
			if (dx < i->first.x) dx = i->first.x;
			if (dy < i->first.y) dy = i->first.y;
			if (dt < i->second.get_deg()) dt = i->second.get_deg();
		}
	}
	if (empty) {
		data.push_back(pol3v1_t());
	}
}
pol3v_t& pol3v_t::operator=(const pol3v_t &x) {
	if (&x != this) {
		data.clear();
		data = x.data;
		dx = x.dx; dy = x.dy; dt = x.dt;
	}
	return *this;
}
pol3v_t pol3v_t::operator+(const pol3v_t &x) const {
	polmap_t tmp;
	for (auto i = data.begin(); i != data.end(); ++i) {
		if (tmp.find(i->first) == tmp.end()) {
			tmp.insert(*i);
		} else {
			tmp[i->first] += i->second;
		}
	}
	for (auto i = x.data.begin(); i != x.data.end(); ++i) {
		if (tmp.find(i->first) == tmp.end()) {
			tmp.insert(*i);
		} else {
			tmp[i->first] += i->second;
		}
	}
	return tmp;
}
pol3v_t pol3v_t::operator-(const pol3v_t &x) const {
	polmap_t tmp;
	for (auto i = data.begin(); i != data.end(); ++i) {
		if (tmp.find(i->first) == tmp.end()) {
			tmp.insert(*i);
		} else {
			tmp[i->first] += i->second;
		}
	}
	for (auto i = x.data.begin(); i != x.data.end(); ++i) {
		if (tmp.find(i->first) == tmp.end()) {
			tmp.insert(pol3v1_t(i->first, pol_t() - i->second));
		} else {
			tmp[i->first] -= i->second;
		}
	}
	return tmp;
}
pol3v_t pol3v_t::operator*(const pol3v_t &x) const {
	polmap_t tmp;
	for (auto i = data.begin(); i != data.end(); ++i) {
		for (auto j = x.data.begin(); j != x.data.end(); ++j) {
			if (tmp.find(i->first + j->first) == tmp.end()) {
				tmp.insert(pol3v1_t(i->first + j->first, i->second * j->second));
			} else {
				tmp[i->first + j->first] += i->second * j->second;
			}
		}
	}
	return tmp;
}
pol3v_t& pol3v_t::operator+=(const pol3v_t &x) {
	return *this = *this + x;
}
pol3v_t& pol3v_t::operator-=(const pol3v_t &x) {
	return *this = *this - x;
}
pol3v_t& pol3v_t::operator*=(const pol3v_t &x) {
	return *this = *this * x;
}
pol_t pol3v_t::subst(const pol_t &x, const pol_t &y) const {
	pol_t res;
	for (auto i = data.begin(), end = data.end(); i != end; ++i) {
		res += x.pow(i->first.x) * y.pow(i->first.y) * i->second;
	}
	return res;
}
pol3v_t pol3v_t::get_same_rand_pol() const {
	pol3v_t res = *this;
	for (auto i = res.data.begin(), end = res.data.end(); i != end; ++i) {
		i->second = get_rand_pol(i->second.get_deg(), false);
	}
	return res;
}
void pol3v_t::read(FILE *fin, pol_format_t format) {
	size_t count = 0;
	data.clear();
	switch (format)
	{
	case BINPR:
	{
		std::fread(&count, 1, sizeof(count), fin);
		if (!count) {
			data.push_back(pol3v1_t());
			dx = 0; dy = 0; dt = 0;
		}
		deg_t deg;
		for (size_t i = 0; i < count; ++i) {
			std::fread(&(deg.key), sizeof(deg.key), 1, fin);
			pol_t tmp;
			tmp.read(fin, BINPR);
			if (tmp) {
				data.push_back(pol3v1_t(deg, tmp));
				if (deg.x > dx) dx = deg.x;
				if (deg.y > dy) dy = deg.y;
				if (tmp.get_deg() > dt) dt = tmp.get_deg();
			}
		}
		break;
	}
	case TEXTPR:
		std::fscanf(fin, "%zu", &count);
		if (!count) {
			data.clear();
			data.push_back(pol3v1_t());
			dx = 0; dy = 0; dt = 0;
		} else {
			for (size_t i = 0; i < count; ++i) {
				u32_t v[2];
				fscanf(fin, "%u%u", v, v + 1);
				pol_t tmp;
				tmp.read(fin, TEXTPR);
				if (tmp) {
					data.push_back(pol3v1_t(deg_t(v[0], v[1]), tmp));
					if (v[0] > dx) dx = v[0];
					if (v[1] > dy) dy = v[1];
					if (tmp.get_deg() > dt) dt = tmp.get_deg();
				}
			}
		}
		break;
	case HTMLPR:
		break;
	default:
		break;
	}
}
void pol3v_t::write(FILE *fout, pol_format_t format, bool eol) const {
	size_t count = data.size();
	switch (format)
	{
	case pol::BINPR:
		std::fwrite(&count, 1, sizeof(count), fout);
		for (auto i = data.begin(), end = data.end(); i != end; ++i) {
			std::fwrite(&(i->first.key), sizeof(i->first.key), 1, fout);
			i->second.write(fout, BINPR);
		}
		break;
	case TEXTPR:
		std::fprintf(fout, "%zu\n", count);
		for (auto i = data.begin(), end = data.end(); i != end; ++i) {
			fprintf(fout, "%u %u ", i->first.x, i->first.y);
			i->second.write(fout, TEXTPR, true);
		}
		if (eol) {
			fprintf(fout, "\n");
		}
		break;
	case pol::HTMLPR:
		if (!count) {
			fprintf(fout, "0");
		} else {
			auto i = data.rbegin();
			if (i->first.x > 1) {
				fprintf(fout, "x<sup>%u</sup>", i->first.x);
			}
			else if (i->first.x == 1) {
				fprintf(fout, "x");
			}
			if (i->first.y > 1) {
				fprintf(fout, "y<sup>%u</sup>", i->first.y);
			}
			else if (i->first.y == 1) {
				fprintf(fout, "y");
			}
			fprintf(fout, "(");
			i->second.write(fout, HTMLPR);
			fprintf(fout, ")");
			++i;
			for (auto end = data.rend(); i != end; ++i) {
				fprintf(fout, " + ");
				if (i->first.x > 1) {
					fprintf(fout, "x<sup>%u</sup>", i->first.x);
				}
				else if (i->first.x == 1) {
					fprintf(fout, "x");
				}
				if (i->first.y > 1) {
					fprintf(fout, "y<sup>%u</sup>", i->first.y);
				}
				else if (i->first.y == 1) {
					fprintf(fout, "y");
				}
				fprintf(fout, "(");
				i->second.write(fout, HTMLPR);
				fprintf(fout, ")");
			}
		}
		if (eol) {
			fprintf(fout, "<br>");
		}
		break;
	default:
		break;
	}
}

File::File() :
	file(NULL),
	opened(false),
	mode(READ),
	buf(0),
	pos(0),
	count(0)
{
}
File::File(FILE *_file, file_mode_t m) :
file(_file),
opened(false),
mode(m),
buf(0),
pos(0),
count(0)
{
}
File::File(const char *file_name, file_mode_t m) :
file(NULL),
opened(false),
mode(m),
buf(0),
pos(0),
count(0)
{
	if (!open(file_name, m)) {
		throw "Can't open the file";
	}
}
File::~File() {
	close();
}
bool File::open(const char *file_name, file_mode_t m) {
	if ((mode = m) == READ) {
		file = fopen(file_name, "rb");
	} else {
		file = fopen(file_name, "wb");
	}
	return opened = !!(file);
}
void File::close() {
	if (opened) {
		flush();
		fclose(file);
		opened = false;
		mode = READ;
	}
}
void File::flush() {
	if (mode == WRITE && count) {
		buf = (buf >> pos) & ((1ull << count) - 1);
		fwrite(&buf, 1, (count - 1) / 8 + 1, file);
	}
	count = 0;
	pos = 0;
	buf = 0;
}
size_t File::read(u61_t &block) {
	if (mode != READ || !opened) {
		throw "Incorrect read";
	}
	u64_t res = 0;
	size_t cur = 0;
	if (count) {
		res = (buf >> pos) & ((1ull << MODSIZE) - 1) & ((1ull << count) - 1);
		cur = count < MODSIZE ? count : MODSIZE;
	}
	if (cur < MODSIZE) {
		buf = 0;
		count = 8 * fread(&buf, 1, (MODSIZE - cur - 1) / 8 + 1, file);
		read_size += count / 8;
		res |= (buf << cur) & ((1ull << MODSIZE) - 1);
	}
	if (res >= MOD || (res < (1ull << (MODSIZE - 1)))) {
		res &= (1ull << (MODSIZE - 1)) - 1;
		if (count > MODSIZE - 1 - cur) {
			count -= MODSIZE - 1 - cur;
			pos = MODSIZE - 1 - cur;
			cur = MODSIZE - 1;
		} else {
			cur += count;
			count = 0;
		}
	} else {
		count -= MODSIZE - cur;
		pos = MODSIZE - cur;
		cur = MODSIZE;
	}
	block = res;
	return cur;
}
size_t File::write(const u61_t &block, size_t block_size) {
	if (mode != WRITE || !opened) {
		throw "Incorrect write";
	}
	u64_t tmp = block;
	size_t size;
	if (tmp >= (1ull << (MODSIZE - 1))) {
		size = MODSIZE;
	} else {
		size = block_size;
	}
	if (!count) {
		buf = 0;
	}
	buf |= tmp << count;
	if (size + count >= 8 * sizeof(buf)) {
		fwrite(&buf, 1, sizeof(buf), file);
	}
	if (size + count > 8 *sizeof(buf)) {
		buf = tmp >> (sizeof(buf)* 8 - count);
	}
	count = (size + count) % (8 * sizeof(buf));
	return size;
}

random_t::random_t() :
	generator(std::time(NULL))
{
}
u61_t random_t::operator()() {
	return u64_t((generator() - 1.) / generator.max() * MOD);
}

pol_t pol::gcd(const pol_t &x, const pol_t &y) {
	pol_t a = x, b = y, q;
	while (q = a % b, q) {
		a = b;
		b = q;
	}
	return b.normalize();
}
pol_t pol::get_x() {
	pol_t res;
	res.change_degree(1);
	res[1] = 1;
	return res;
}

pol_t pol::get_rand_pol(u32_t deg, bool flag) {
	static random_t get_coef;
	pol_t res;
	res.change_degree(deg);
	//res[deg] = flag ? 1 : get_coef();
	while (!res(deg)) {
		res[deg] = get_coef();
	}
	for (u32_t i = 0; i < deg; ++i) {
		res[i] = get_coef();
	}
	return res;
}

pol3v_t pol::XXX_get_f(const u64_t max_deg, const u32_t key_deg) {
	polmap_t tmp;
	for (u32_t i = 0; i <= max_deg; ++i) {
		tmp.insert(pol3v1_t(deg_t(i, i), get_rand_pol(key_deg - 1)));
	}
	return tmp;
}
pol3v_t pol::XXX_get_mes(const u64_t max_deg, const u32_t key_deg, const u64_t size, const std::vector<u61_t> &m) {
	polmap_t tmp;
	pol_t tmp2;
	u64_t block_count = m.size();
	tmp2.change_degree(key_deg - 1);
	for (u32_t i = 0; i < max_deg; ++i) {
		for (u32_t j = 0; j < key_deg; ++j) {
			if (i * key_deg + j < block_count) {
				tmp2[j] = m[block_count - 1 - (i * key_deg + j)];
			}
			else {
				tmp2(j) = 0;
			}
		}
		tmp.insert(pol3v1_t(deg_t(i, i), tmp2));
	}
	tmp2.change_degree(0);
	tmp2[0] = size;
	tmp.insert(pol3v1_t(deg_t(u32_t(max_deg), u32_t(max_deg)), tmp2));
	return tmp;
}
void pol::XXX_get_fact(const vector_t &fact, u64_t deg, vector_t &res) {
	size_t num = fact.size();
	std::vector<std::list<pol_t>> tmp;
	tmp.resize(deg + 1);
	for (size_t i = 0; i < num; ++i) {
		u32_t pol_deg = fact[i].get_deg();
		for (auto j = tmp.rbegin(); j != tmp.rend(); ++j) {
			for (auto it = j->begin(); it != j->end(); ++it) {
				if (pol_deg + it->get_deg() <= deg) {
					tmp[pol_deg + it->get_deg()].push_back(fact[i] * *it);
				}
			}
		}
		if (pol_deg <= deg) {
			tmp[pol_deg].push_back(fact[i]);
		}
	}
	res.clear();
	for (auto i = tmp[deg].begin(); i != tmp[deg].end(); ++i) {
		res.push_back(*i);
	}
}
