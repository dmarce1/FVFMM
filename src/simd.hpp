/*
 * simd_vector.hpp
 *
 *  Created on: Jun 7, 2015
 *      Author: dmarce1
 */

#ifndef SIMD_VECTOR_HPP_
#define SIMD_VECTOR_HPP_

#include <cstdlib>
#include "immintrin.h"

const std::size_t simd_len = 4;

class simd_vector {
private:
	__m256d v;
public:
	inline simd_vector() = default;
	inline ~simd_vector() = default;
	inline simd_vector(const simd_vector&) = default;
	inline simd_vector(const __m256d& mm) {
		v = mm;
	}
	inline simd_vector(double d) {
		v =_mm256_set_pd(d,d,d,d);
	}
	inline double sum() const {
		double r = ZERO;
		for( integer i = 0; i != simd_len; ++i) {
			r += (*this)[i];
		}
		return r;
	}
	inline simd_vector& operator=(const __m256d& mm) {
		v = mm;
		return *this;
	}
	inline simd_vector(simd_vector&&) = default;
	inline simd_vector& operator=(const simd_vector&) = default;
	inline simd_vector& operator=(simd_vector&&) = default;

	inline simd_vector operator+( const simd_vector& other) const {
		return simd_vector(_mm256_add_pd(v, other.v));
	}
	inline simd_vector operator-( const simd_vector& other) const {
		return simd_vector(_mm256_sub_pd(v, other.v));
	}
	inline simd_vector operator*( const simd_vector& other) const {
		return simd_vector(_mm256_mul_pd(v, other.v));
	}
	inline simd_vector operator/( const simd_vector& other) const {
		return simd_vector(_mm256_div_pd(v, other.v));
	}
	inline simd_vector operator+() const {
		return *this;
	}
	inline simd_vector operator-() const {
		return simd_vector(ZERO) - *this;
	}
	inline simd_vector& operator+=( const simd_vector& other) {
		*this = *this + other;
		return *this;
	}
	inline simd_vector& operator-=( const simd_vector& other) {
		*this = *this - other;
		return *this;
	}
	inline simd_vector& operator*=( const simd_vector& other) {
		*this = *this * other;
		return *this;
	}
	inline simd_vector& operator/=( const simd_vector& other) {
		*this = *this / other;
		return *this;
	}

	inline simd_vector operator*( double d) const {
		const simd_vector other = d;
		return other * *this;
	}
	inline simd_vector operator/(double d) const {
		const simd_vector other = d;
		return *this / other;
	}

	inline simd_vector operator*=( double d) {
		*this = *this * d;
		return *this;
	}
	inline simd_vector operator/=(double d) {
		*this = *this / d;
		return *this;
	}
	inline double& operator[](std::size_t i) {
		double* a = reinterpret_cast<double*>(&v);
		return a[i];
	}
	inline double operator[](std::size_t i) const {
		const double* a = reinterpret_cast<const double*>(&v);
		return a[i];
	}

	friend simd_vector sqrt(const simd_vector&);
	friend simd_vector operator*(double, const simd_vector& other);
	friend simd_vector operator/(double, const simd_vector& other);
};

inline simd_vector sqrt(const simd_vector& v) {
	return simd_vector(_mm256_sqrt_pd(v.v));
}

inline simd_vector operator*(double d, const simd_vector& other) {
	const simd_vector a = d;
	return a * other;
}

inline simd_vector operator/(double d, const simd_vector& other) {
	const simd_vector a = d;
	return a / other;
}

#endif /* SIMD_VECTOR_HPP_ */
