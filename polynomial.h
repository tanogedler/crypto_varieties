// polynomial.h

#pragma once

#include <cstdio>
#include <utility>
#include <vector>
#include <map>
#include <list>
#include <random>
#include <ctime>
#include "hash.h"

namespace pol {

	typedef unsigned long long u64_t;
	typedef unsigned u32_t;
	typedef unsigned char u8_t;

	const u32_t POL1V_MAXLENGTH = 4096;
	const u32_t M32 = 0xffffffff;
	const u32_t KEY_DEG = 3;
	const u32_t MODSIZE = 61;
	const u32_t BUFSIZE = 8;
	const u64_t MOD = 2305843009213693951; // 2147483647; //  // 2^61 - 1

	class u61_t {
	private:
		u64_t data;
	public:
		u61_t(u64_t x = 0);
		 operator u64_t() const;
		 u61_t operator+(const u61_t &) const;
	    u61_t operator-(const u61_t &) const;
		u61_t operator*(const u61_t &) const;
		u61_t operator/(const u61_t &) const;
		u61_t& operator+=(const u61_t &);
		u61_t& operator-=(const u61_t &);
		u61_t& operator*=(const u61_t &);
		u61_t& operator/=(const u61_t &);
		u61_t inverse() const;
		u61_t pow(u32_t) const;
		inline u64_t* get_addr();
		inline const u64_t* get_addr() const;
		inline u64_t& get_buf();
		inline const u64_t& get_buf() const;
		void print(FILE *, u32_t) const;
	};

	enum pol_format_t {
		BINPR,
		TEXTPR,
		HTMLPR
	};

	class pol_t {
	private:
		std::vector<u61_t> data;
		u32_t degree;
	public:
		typedef std::vector<pol_t> vector_t;
		pol_t();
		pol_t& operator=(const pol_t&);
		void change_degree(u32_t);
		pol_t& operator+=(const pol_t&);
		pol_t& operator-=(const pol_t&);
		pol_t& operator*=(const pol_t&);
		pol_t& operator%=(const pol_t&);
		pol_t& operator/=(const pol_t&);
		pol_t operator+(const pol_t&) const;
		pol_t operator-(const pol_t&) const;
		pol_t operator*(const pol_t&) const;
		pol_t operator%(const pol_t&) const;
		pol_t operator/(const pol_t&) const;
		operator bool() const;
		bool operator==(const pol_t &) const;
		pol_t pow(u32_t) const;
		pol_t mpow(u64_t, const pol_t &) const;
		pol_t mpow(u64_t, u64_t, const pol_t &) const;
		pol_t normalize() const;
		pol_t deriv() const;
		void read(FILE * fin, pol_format_t format = TEXTPR);
		void write(FILE * fout, pol_format_t format = HTMLPR, bool endofline = true) const;
		u61_t& operator[](u32_t);
		const u61_t& operator[](u32_t) const;
		u64_t& operator()(u32_t);
		const u64_t& operator()(u32_t) const;
		 u32_t get_deg() const;
		void fact(u32_t, vector_t &) const;
		void ddf(u32_t, vector_t &) const;
		void split(u32_t, vector_t &) const;
	};

	union deg_t {
		u64_t key;
		struct {
			u32_t x, y;
		};
		deg_t();
		deg_t(u32_t, u32_t);
		bool operator<(const deg_t &) const;
		bool operator==(const deg_t &) const;
		deg_t operator+(const deg_t&) const;
	};

	typedef std::pair<deg_t, pol_t> pol3v1_t;
	typedef std::list<pol3v1_t> list_t;
	typedef std::map<deg_t, pol_t> polmap_t;
	typedef pol_t::vector_t vector_t;

	class pol3v_t {
	private:
		list_t data;
		u32_t dx;
		u32_t dy;
		u32_t dt;
	public:
		pol3v_t();
		pol3v_t(const polmap_t &);
		pol3v_t& operator=(const pol3v_t &);
		pol3v_t operator+(const pol3v_t&) const;
		pol3v_t operator-(const pol3v_t&) const;
		pol3v_t operator*(const pol3v_t&) const;
		pol3v_t& operator+=(const pol3v_t &);
		pol3v_t& operator-=(const pol3v_t &);
		pol3v_t& operator*=(const pol3v_t &);
		pol_t subst(const pol_t &, const pol_t &) const;
		pol3v_t get_same_rand_pol() const;
		void read(FILE *, pol_format_t);
		void write(FILE * fout, pol_format_t format = HTMLPR, bool endofline = 1) const;
	};

	enum file_mode_t {
		READ, WRITE
	};

	class File {
	private:
		FILE *file;
		bool opened;
		file_mode_t mode;
		u64_t buf;
		size_t pos, count, read_size;
	public:
		File();
		File(FILE *, file_mode_t);
		File(const char *, file_mode_t);
		~File();
		bool open(const char *, file_mode_t);
		void close();
		void flush();
		size_t read(u61_t &);
		size_t write(const u61_t &, size_t = MODSIZE - 1);
	};

	class random_t {
	private:
		std::mt19937_64 generator;
	public:
		random_t();
		u61_t operator()();
	};

	pol_t gcd(const pol_t &, const pol_t &);
	pol_t get_x();
	pol_t get_rand_pol(u32_t deg, bool f = true);
	template <unsigned CLOSE_KEY_DEG = KEY_DEG, unsigned OPEN_KEY_DEG = KEY_DEG>
	pol3v_t get_open_key(const pol_t &, const pol_t &);
	template <unsigned Key_Deg = 2 * KEY_DEG>
	void encrypt(File &, FILE *, const pol3v_t &);
	template <unsigned CLOSE_KEY_DEG = KEY_DEG, unsigned OPEN_KEY_DEG = KEY_DEG>
	void decrypt(FILE *, File &, const pol_t &, const pol_t &);

	pol3v_t XXX_get_f(const u64_t max_deg, const u32_t Key_Deg);
	pol3v_t XXX_get_mes(const u64_t max_deg, const u32_t Key_Deg, const u64_t size, const std::vector<u61_t> &);
	void XXX_get_fact(const vector_t &fact, u64_t deg, vector_t &res);

	#include "polinomial.hpp"
};
