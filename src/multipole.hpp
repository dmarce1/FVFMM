#ifndef multipole_H_
#define multipole_H_

#include "defs.hpp"
#include "space_vector.hpp"

class multipole: public std::array<real,MP> {
private:
	const real delta[3][3] = { { real(1.0), real(0.0), real(0.0) }, { real(0.0), real(1.0), real(0.0) }, { real(0.0),
			real(0.0), real(1.0) } };
	const size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	const size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } },
			{ { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8, 9 } } };
public:
	multipole();
	real operator ()() const;
	real& operator ()();
	real operator ()(integer i) const;
	real& operator ()(integer i);
	real operator ()(integer i, integer j) const;
	real& operator ()(integer i, integer j);
	real operator ()(integer i, integer j, integer k) const;
	real& operator ()(integer i, integer j, integer k);
	multipole& operator =(const multipole& expansion);
	multipole& operator =(real expansion);
	multipole operator>>(const space_vector& dX) const;
	multipole& operator>>=(const space_vector& Y);
	std::array<real, MP>& operator +=(const std::array<real, MP>& vec);
	~multipole();
};

#endif /* multipole_H_ */
