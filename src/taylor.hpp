/*
 * taylor.hpp
 *
 *  Created on: Jun 9, 2015
 *      Author: dmarce1
 */

#ifndef TAYLOR_HPP_
#define TAYLOR_HPP_

#include "defs.hpp"

#include <array>
#include <cmath>

#define MAX_ORDER 5

struct taylor_consts {
	static const real delta[3][3];
	static integer map2[3][3];
	static integer map3[3][3][3];
	static integer map4[3][3][3][3];
};

template<int npoles, class T = real>
class taylor {
private:
	static constexpr integer sizes[MAX_ORDER] = {1, 4, 10, 20, 35}; //
	static constexpr integer size = sizes[npoles-1];
	static taylor_consts tc;
	std::array<T, size> data;
public:
	taylor() = default;
	~taylor() = default;
	taylor(const taylor&) = default;
	taylor(taylor&&) = default;
	taylor& operator=(const taylor&) = default;
	taylor& operator=(taylor&&) = default;

	taylor& operator=(T d) {
		for( integer i = 0; i != size; ++i ) {
			data[i] = d;
		}
		return *this;
	}

	taylor& operator*=(T d) {
		for( integer i = 0; i != size; ++i ) {
			data[i] *= d;
		}
		return *this;
	}

	taylor& operator/=(T d) {
		for( integer i = 0; i != size; ++i ) {
			data[i] /= d;
		}
		return *this;
	}

	taylor& operator+=(const taylor& other) {
		for( integer i = 0; i != size; ++i) {
			data[i] += other.data[i];
		}
		return *this;
	}

	taylor& operator-=(const taylor& other) {
		for( integer i = 0; i != size; ++i) {
			data[i] -= other.data[i];
		}
		return *this;
	}

	taylor operator+(const taylor& other) const {
		taylor r = *this;
		r += other;
		return r;
	}

	taylor operator-(const taylor& other) const {
		taylor r = *this;
		r -= other;
		return r;
	}

	taylor operator*(const T& d) const {
		taylor r = *this;
		r *= d;
		return r;
	}

	taylor operator/(const T& d) const {
		taylor r = *this;
		r /= d;
		return r;
	}

	taylor operator+() const {
		return *this;
	}

	taylor operator-() const {
		taylor r = *this;
		for( integer i = 0; i != size; ++i) {
			r.data[i] = -r.data[i];
		}
		return r;
	}

	T operator()() const {
		return data[0];
	}

	T operator()(integer i) const {
		return data[1 + i];
	}

	T operator()(integer i, integer j) const {
		return data[tc.map2[i][j]];
	}

	T operator()(integer i, integer j, integer k) const {
		return data[tc.map3[i][j][k]];
	}

	T operator()(integer i, integer j, integer k, integer l) const {
		return data[tc.map4[i][j][k][l]];
	}

	T& operator()() {
		return data[0];
	}

	T& operator()(integer i) {
		return data[1 + i];
	}

	T& operator()(integer i, integer j) {
		return data[tc.map2[i][j]];
	}

	T& operator()(integer i, integer j, integer k) {
		return data[tc.map3[i][j][k]];
	}

	T& operator()(integer i, integer j, integer k, integer l) {
		return data[tc.map4[i][j][k][l]];
	}

	taylor& operator>>=(const std::array<T,NDIM>& X) {
		const taylor& A = *this;
		taylor B = A;

		if( npoles > 1 ) {
			for( integer a = 0; a < NDIM; ++a) {
				B(a) += A() * X[a];
			}
			if( npoles > 2 ) {
				for( integer a = 0; a < NDIM; ++a) {
					for( integer b = a; b < NDIM; ++b) {
						B(a,b) += A(a) * X[b] + X[a] * A(b);
						B(a,b) += A()*X[a]*X[b];
					}
				}
				if( npoles > 3 ) {
					for( integer a = 0; a < NDIM; ++a) {
						for( integer b = a; b < NDIM; ++b) {
							for( integer c = b; c < NDIM; ++c) {
								B(a,b,c) += A(a,b) * X[c] + A(b,c) * X[a] + A(c,a) * X[b];
								B(a,b,c) += A(a) * X[b] * X[c] + A(b) * X[c] * X[a] + A(c) * X[a] * X[b];
								B(a,b,c) += A() * X[a] * X[b] * X[c];
							}
						}
					}
				}
			}
		}
		*this = B;
		return *this;
	}

	taylor operator>>(const std::array<T,NDIM>& X) const {
		taylor r = *this;
		r >>= X;
		return r;
	}

	taylor& operator<<=(const std::array<T,NDIM>& X ) {
		const taylor& A = *this;
		taylor B = A;

		if( npoles > 1 ) {
			for (integer a = 0; a != NDIM; a++) {
				B() += A(a) * X[a];
			}
			if( npoles > 2 ) {
				for (integer a = 0; a != NDIM; a++) {
					for (integer b = 0; b != NDIM; b++) {
						B() += A(a, b) * X[a] * X[b] / T(2);
					}
				}
				for (integer a = 0; a != NDIM; a++) {
					for (integer b = 0; b != NDIM; b++) {
						B(a) += A(a, b) * X[b];
					}
				}
				if( npoles > 3 ) {
					for (integer a = 0; a != NDIM; a++) {
						for (integer b = 0; b != NDIM; b++) {
							for (integer c = 0; c != NDIM; c++) {
								B() += A(a, b, c) * X[a] * X[b] * X[c] / T(6);
							}
						}
					}
					for (integer a = 0; a != NDIM; a++) {
						for (integer b = 0; b != NDIM; b++) {
							for (integer c = 0; c != NDIM; c++) {
								B(a) += A(a, b, c) * X[b] * X[c] / T(2);
							}
						}
					}
					for (integer a = 0; a != NDIM; a++) {
						for (integer b = 0; b != NDIM; b++) {
							for (integer c = a; c != NDIM; c++) {
								B(a, c) += A(a, b, c) * X[b];
							}
						}
					}
				}
			}
		}
		*this = B;
		return *this;
	}

	taylor operator<<(const std::array<T,NDIM>& X) const {
		taylor r = *this;
		r <<= X;
		return r;
	}

	void set_basis(const std::array<T,NDIM>& X) {

		taylor& A = *this;
		T r = ZERO;
		for (integer d = 0; d != NDIM; ++d) {
			r += X[d] * X[d];
		}
		const T r2inv = T(1) / r;

		const T d0 = -sqrt(r2inv);
		A() = d0;

		if( npoles > 1) {
			const T d1 = -d0 * r2inv;
			for (integer a = 0; a != NDIM; a++) {
				A(a) = X[a] * d1;
			}
			if( npoles > 2 ) {
				const T d2 = -T(3) * d1 * r2inv;
				for (integer a = 0; a != NDIM; a++) {
					for (integer b = a; b != NDIM; b++) {
						A(a, b) = X[a] * X[b] * d2;
						A(a, b) += tc.delta[a][b] * d1;
					}
				}
				if( npoles > 3 ) {
					const T d3 = -T(5) * d2 * r2inv;
					for (integer a = 0; a != NDIM; a++) {
						for (integer b = a; b != NDIM; b++) {
							for (integer c = b; c != NDIM && b != NDIM; c++) {
								A(a, b, c) = X[a] * X[b] * X[c] * d3;
								A(a, b, c) += (tc.delta[a][b] * X[c] + tc.delta[b][c] * X[a] + tc.delta[c][a] * X[b]) * d2;
							}
						}
					}
					if( npoles > 4 ) {
						const T d4 = -T(7) * d3 * r2inv;
						for (integer a = 0; a != NDIM; a++) {
							for (integer b = a; b != NDIM; b++) {
								for (integer c = b; c != NDIM; c++) {
									for( integer d = c; d != NDIM && c != NDIM; ++d) {
										A(a,b,c,d) = X[a] * X[b] * X[c] * X[d] * d4;
										A(a,b,c,d) += (tc.delta[a][b] * X[c] * X[d]
												+ tc.delta[a][c] * X[b] * X[d]
												+ tc.delta[a][d] * X[b] * X[c]
												+ tc.delta[b][c] * X[a] * X[d]
												+ tc.delta[b][d] * X[a] * X[c]
												+ tc.delta[c][d] * X[a] * X[b]) * d3;
										A(a,b,c,d) += (tc.delta[a][b] * tc.delta[c][d]
												+ tc.delta[a][c] * tc.delta[b][d]
												+ tc.delta[a][d] * tc.delta[b][c]) * d2;

									}
								}
							}
						}
					}
				}
			}
		}
	}

};

template<int npoles, class T>
taylor_consts taylor<npoles, T>::tc;

using multipole = taylor<4,real>;
using expansion = taylor<4,real>;
#endif /* TAYLOR_HPP_ */
