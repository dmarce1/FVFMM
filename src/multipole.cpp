#include "multipole.hpp"

multipole::multipole() {
}

real multipole::operator ()() const {
	return (*this)[0];
}
real& multipole::operator ()() {
	return (*this)[0];
}

real multipole::operator ()(integer i) const {
	return (*this)[1 + i];
}

real& multipole::operator ()(integer i) {
	return (*this)[1 + i];
}

real multipole::operator ()(integer i, integer j) const {
	return (*this)[4 + map2[i][j]];
}

real& multipole::operator ()(integer i, integer j) {
	return (*this)[4 + map2[i][j]];
}

real multipole::operator ()(integer i, integer j, integer k) const {
	return (*this)[10 + map3[i][j][k]];
}

real& multipole::operator ()(integer i, integer j, integer k) {
	return (*this)[10 + map3[i][j][k]];
}

multipole& multipole::operator =(const multipole& expansion) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] = expansion[i];
	}
	return *this;
}

multipole& multipole::operator =(real expansion) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] = expansion;
	}
	return *this;
}

multipole multipole::operator>>(const space_vector& dX) const {
	multipole you = *this;
	you >>= dX;
	return you;
}

multipole& multipole::operator>>=(const space_vector& Y) {
	multipole& me = *this;
	multipole tmp;
	tmp = me;
	for (int p = 0; p < 3; p++) {
		me(p) += tmp() * Y[p];
	}
	for (int p = 0; p < 3; p++) {
		for (int q = p; q < 3; q++) {
			me(p, q) += tmp(q) * Y[p] + tmp(p) * Y[q];
			me(p, q) += tmp() * Y[p] * Y[q];
		}
	}
	for (int p = 0; p < 3; p++) {
		for (int q = p; q < 3; q++) {
			for (int r = q; r < 3; r++) {
				me(p, q, r) += tmp(p) * Y[q] * Y[r] + tmp(q) * Y[r] * Y[p] + tmp(r) * Y[p] * Y[q];
				me(p, q, r) += tmp() * Y[p] * Y[q] * Y[r] + tmp(p, q) * Y[r] + tmp(q, r) * Y[p] + tmp(r, p) * Y[q];
			}
		}
	}
	return me;
}

std::array<real, MP>& multipole::operator +=(const std::array<real, MP>& vec) {
	for (integer i = 0; i < MP; i++) {
		(*this)[i] += vec[i];
	}
	return *this;
}

multipole::~multipole() {
}
