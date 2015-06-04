#include "expansion.hpp"
#include <cmath>

expansion::expansion() {
}

real expansion::operator ()() const {
	return (*this)[0];
}
real& expansion::operator ()() {
	return (*this)[0];
}

real expansion::operator ()(int i) const {
	return (*this)[1 + i];
}
real& expansion::operator ()(int i) {
	return (*this)[1 + i];
}

real expansion::operator ()(int i, int j) const {
	return (*this)[4 + map2[i][j]];
}
real& expansion::operator ()(int i, int j) {
	return (*this)[4 + map2[i][j]];
}

real expansion::operator ()(int i, int j, int k) const {
	return (*this)[10 + map3[i][j][k]];
}
real& expansion::operator ()(int i, int j, int k) {
	return (*this)[10 + map3[i][j][k]];
}

expansion& expansion::operator =(const expansion& expansion) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] = expansion[i];
	}
	return *this;
}

expansion& expansion::operator =(real expansion) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] = expansion;
	}
	return *this;
}

expansion expansion::operator<<(const space_vector& dX) const {
	expansion you = *this;
	you <<= dX;
	return you;
}

expansion& expansion::operator<<=(const space_vector& dX) {
	expansion& me = *this;
	for (integer a = 0; a < 3; a++) {
		me() += me(a) * dX[a];
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			me() += me(a, b) * dX[a] * dX[b] * 0.5;
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			me(a) += me(a, b) * dX[b];
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				me() += me(a, b, c) * dX[a] * dX[b] * dX[c] * (1.0 / 6.0);
			}
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				me(a) += me(a, b, c) * dX[b] * dX[c] * 0.5;
			}
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = a; c < 3; c++) {
				me(a, c) += me(a, b, c) * dX[b];
			}
		}
	}
	return me;
}

real expansion::translate_to_particle(const space_vector& dX) const {
	const auto& L = *this;
	real this_phi = L();
	for (integer a = 0; a < 3; a++) {
		this_phi += L(a) * dX[a];
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			this_phi += L(a, b) * dX[a] * dX[b] * 0.5;
		}
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = 0; b < 3; b++) {
			for (integer c = 0; c < 3; c++) {
				this_phi += L(a, b, c) * dX[a] * dX[b] * dX[c] * (1.0 / 6.0);
			}
		}
	}
	return this_phi;
}

std::array<real, LP>& expansion::operator +=(const std::array<real, LP>& vec) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] += vec[i];
	}
	return *this;
}

std::array<real, LP>& expansion::operator -=(const std::array<real, LP>& vec) {
	for (int i = 0; i < LP; i++) {
		(*this)[i] -= vec[i];
	}
	return *this;
}

void expansion::compute_D(const space_vector& Y) {
	expansion& me = *this;
	real y0 = 0.0;
	for (integer d = 0; d != NDIM; ++d) {
		y0 += Y[d] * Y[d];
	}
	const real r2inv = 1.0 / y0;
	const real d0 = -std::sqrt(r2inv);
	const real d1 = -d0 * r2inv;
	const real d2 = -3.0 * d1 * r2inv;
	const real d3 = -5.0 * d2 * r2inv;
	me() = d0;
	for (integer a = 0; a < 3; a++) {
		me(a) = Y[a] * d1;
		for (integer b = a; b < 3; b++) {
			me(a, b) = Y[a] * Y[b] * d2 + delta[a][b] * d1;
			for (integer c = b; c < 3; c++) {
				me(a, b, c) = Y[a] * Y[b] * Y[c] * d3;
				me(a, b, c) += (delta[a][b] * Y[c] + delta[b][c] * Y[a] + delta[c][a] * Y[b]) * d2;
			}
		}
	}
}

void expansion::invert() {
	expansion& me = *this;
	for (integer a = 0; a < 3; a++) {
		me(a) = -me(a);
	}
	for (integer a = 0; a < 3; a++) {
		for (integer b = a; b < 3; b++) {
			for (integer c = b; c < 3; c++) {
				me(a, b, c) = -me(a, b, c);
			}
		}
	}
}

expansion::~expansion() {
}

expansion operator*(const multipole& M, const expansion& D) {
	expansion me;
	me() = M() * D();
	for (int a = 0; a < 3; a++) {
		me() -= M(a) * D(a);
	}
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			me() += M(a, b) * D(a, b) * 0.5;
		}
	}
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			for (int c = 0; c < 3; c++) {
				me() -= M(a, b, c) * D(a, b, c) * (1.0 / 6.0);
			}
		}
	}

	for (int a = 0; a < 3; a++) {
		me(a) = M() * D(a);
	}
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			me(a) -= M(a) * D(a, b);
		}
	}
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			for (int c = 0; c < 3; c++) {
				me(a) += M(c, b) * D(a, b, c) * 0.5;
			}
		}
	}

	for (int a = 0; a < 3; a++) {
		for (int b = a; b < 3; b++) {
			me(a, b) = M() * D(a, b);
		}
	}
	for (int a = 0; a < 3; a++) {
		for (int b = a; b < 3; b++) {
			for (int c = 0; c < 3; c++) {
				me(a, b) -= M(c) * D(a, b, c);
			}
		}
	}

	for (int a = 0; a < 3; a++) {
		for (int b = a; b < 3; b++) {
			for (int c = b; c < 3; c++) {
				me(a, b, c) = M() * D(a, b, c);
			}
		}
	}
	return me;
}

std::array<expansion, NDIM> expansion::get_derivatives(const multipole& M, const multipole& N,
		const space_vector& R) const {
	std::array<expansion, NDIM> D;
	for (integer i = 0; i != NDIM; ++i) {
		D[i]() = (*this)(i);
		for (integer a = 0; a < NDIM; ++a) {
			D[i](a) = (*this)(i, a);
			for (integer b = a; b < NDIM; ++b) {
				D[i](a, b) = (*this)(i, a, b);
				for (integer c = b; c != NDIM; ++c) {
					D[i](a, b, c) = real(0.0);
				}
			}
		}
	}
	real r = 0.0;
	for (integer d = 0; d != NDIM; ++d) {
		r += R[d] * R[d];
	}
	r = std::sqrt(r);
	const real rinv = real(1) / r;
	const real D2 = -real(3.0) * std::pow(rinv, real(5));
	const real D3 = +real(15.) * std::pow(rinv, real(7));
	const real D4 = -real(105) * std::pow(rinv, real(9));
	std::array<expansion, NDIM> E1;
	std::array<expansion, NDIM> E2;
	for (integer i = 0; i != NDIM; ++i) {
		E1[i] = 0.0;
		E2[i] = 0.0;
		for (integer j = 0; j != NDIM; ++j) {
			for (integer k = j; k != NDIM; ++k) {
				for (integer l = k; l != NDIM; ++l) {
					E1[i](j, k, l) += delta[i][j] * delta[k][l] * D2;
					E1[i](j, k, l) += delta[i][l] * delta[j][k] * D2;
					E1[i](j, k, l) += delta[i][k] * delta[l][j] * D2;

					E1[i](j, k, l) += delta[i][j] * R[k] * R[l] * D3;
					E1[i](j, k, l) += delta[i][l] * R[j] * R[k] * D3;
					E1[i](j, k, l) += delta[i][k] * R[l] * R[j] * D3;

					E2[i](j, k, l) += R[i] * R[j] * delta[k][l] * D3;
					E2[i](j, k, l) += R[i] * R[l] * delta[j][k] * D3;
					E2[i](j, k, l) += R[i] * R[k] * delta[l][j] * D3;

					E2[i](j, k, l) += R[i] * R[j] * R[k] * R[l] * D4;
				}
			}
		}
	}
	for (integer i = 0; i != NDIM; ++i) {
		for (integer j = 0; j != NDIM; ++j) {
			for (integer k = 0; k != NDIM; ++k) {
				for (integer l = 0; l != NDIM; ++l) {
					D[i]() -= E1[i](j, k, l) * M(j, k, l) * (real(1) / real(6));
					D[i]() += E1[i](j, k, l) * N(j, k, l) * M() / N() * (real(1) / real(6));
					D[i]() -= E2[i](j, k, l) * M(j, k, l) * (real(1) / real(6));
					D[i]() += E2[i](j, k, l) * N(j, k, l) * M() / N() * (real(1) / real(6));
				}
			}
		}
	}

	return D;
}

