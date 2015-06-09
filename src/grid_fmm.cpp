/*
 * grid_fmm.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: dmarce1
 */
#include "grid.hpp"
#include "simd.hpp"
#include "stop_watch.hpp"

void grid::solve_gravity(gsolve_type type) {
	compute_multipoles(type);
	compute_interactions(type);
	compute_expansions(type);
	if (type == RHO) {
		for (integer iii = 0; iii != HN3; ++iii) {
			G[phi_i][iii] = L[0][iii]();
			for (integer d = 0; d != NDIM; ++d) {
				G[gx_i + d][iii] = -L[0][iii](d) - L_c[0][iii](d);
			}
			U[pot_i][iii] = G[phi_i][iii] * U[rho_i][iii];
		}
	} else {
		for (integer iii = 0; iii != HN3; ++iii) {
			dphi_dt[iii] = L[0][iii]();
		}
	}
}

void grid::compute_interactions(gsolve_type type) {

	npair np;
	dpair dp;
	for (integer lev = 0; lev != nlevel; ++lev) {
		std::fill(L[lev].begin(), L[lev].end(), ZERO);
		std::fill(L_c[lev].begin(), L_c[lev].end(), ZERO);
	}

	for (auto iter = ilist_n.begin(); iter != ilist_n.end(); ++iter) {
		std::array < simd_vector, NDIM > X;
		std::array<simd_vector, NDIM> Y;
		taylor<4, simd_vector> m0;
		taylor<4, simd_vector> n0;
		taylor<4, simd_vector> m1;
		taylor<4, simd_vector> n1;
		decltype(iter) oiter = iter;
		for (integer i = 0; i != simd_len && iter != ilist_n.end(); ++i) {
			const integer iii0 = iter->loc.first;
			const integer iii1 = iter->loc.second;
			const integer lev = iter->lev;
			for (integer d = 0; d != NDIM; ++d) {
				X[d][i] = com[lev][iii0][d];
				Y[d][i] = com[lev][iii1][d];
			}
			for (integer j = 0; j != 20; ++j) {
				m0.ptr()[j][i] = M[lev][iii1].ptr()[j];
				if (type == RHO) {
					n0.ptr()[j][i] = M[lev][iii1].ptr()[j] - M[lev][iii0].ptr()[j] * (M[lev][iii1]() / M[lev][iii0]());
				} else {
					n0.ptr()[j][i] = ZERO;
				}
				m1.ptr()[j][i] = M[lev][iii0].ptr()[j];
				if (type == RHO) {
					n1.ptr()[j][i] = M[lev][iii0].ptr()[j] - M[lev][iii1].ptr()[j] * (M[lev][iii0]() / M[lev][iii1]());
				} else {
					n1.ptr()[j][i] = ZERO;
				}
			}
			if (i != simd_len - 1) {
				++iter;
			}
		}
		std::array<simd_vector, NDIM> dX;
		for (integer d = 0; d != NDIM; ++d) {
			dX[d] = X[d] - Y[d];
		}
		taylor<5, simd_vector> D;
		taylor<4, simd_vector> A0, B0, A1, B1;
		D.set_basis(dX);
		B0 = ZERO;
		B1 = ZERO;
		A0() = m0() * D();
		A1() = m1() * D();
		for (integer a = 0; a != NDIM; ++a) {
			A0() -= m0(a) * D(a);
			A1() += m1(a) * D(a);
			for (integer b = 0; b != NDIM; ++b) {
				A0() += m0(a, b) * D(a, b) / real(2);
				A1() += m1(a, b) * D(a, b) / real(2);
				for (integer c = 0; c != NDIM; ++c) {
					A0() -= m0(a, b, c) * D(a, b, c) / real(6);
					A1() += m0(a, b, c) * D(a, b, c) / real(6);
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			A0(a) = +m0() * D(a);
			A1(a) = -m1() * D(a);
			for (integer b = 0; b != NDIM; ++b) {
				A0(a) -= m0(a) * D(a, b);
				A1(a) -= m1(a) * D(a, b);
				for (integer c = 0; c != NDIM; ++c) {
					A0(a) += m0(c, b) * D(a, b, c) / real(2);
					A1(a) -= m1(c, b) * D(a, b, c) / real(2);
					for (integer d = 0; d != NDIM; ++d) {
						B0(a) -= n0(b, c, d) * D(a, b, c, d) / real(6);
						B1(a) -= n1(b, c, d) * D(a, b, c, d) / real(6);
					}
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			for (integer b = a; b != NDIM; ++b) {
				A0(a, b) = m0() * D(a, b);
				A1(a, b) = m1() * D(a, b);
				for (integer c = 0; c != NDIM; ++c) {
					A0(a, b) -= m0(c) * D(a, b, c);
					A1(a, b) += m1(c) * D(a, b, c);
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			for (integer b = a; b != NDIM; ++b) {
				for (integer c = b; c != NDIM; ++c) {
					A0(a, b, c) = +m0() * D(a, b, c);
					A1(a, b, c) = -m1() * D(a, b, c);
				}
			}
		}
		for (integer i = 0; i != simd_len && oiter != ilist_n.end(); ++i, ++oiter) {
			const integer iii0 = oiter->loc.first;
			const integer iii1 = oiter->loc.second;
			const integer lev = oiter->lev;
			for (integer j = 0; j != 20; ++j) {
				L[lev][iii0].ptr()[j] += A0.ptr()[j][i];
				L[lev][iii1].ptr()[j] += A1.ptr()[j][i];
				if (type == RHO) {
					L_c[lev][iii0].ptr()[j] += B0.ptr()[j][i];
					L_c[lev][iii1].ptr()[j] += B1.ptr()[j][i];
				}
			}
		}
	}
	for (auto iter = ilist_d.begin(); iter != ilist_d.end(); ++iter) {
		const integer iii0 = iter->first;
		const integer iii1 = iter->second;
		const integer lev = 0;
		space_vector dX = (com[lev][iii0] - com[lev][iii1]);
		const real r = dX.abs();
		const real rinv = ONE / r;
		const real r3inv = ONE / (r * r * r);
		L[0][iii0]() -= M[0][iii1]() * rinv;
		if (type == RHO) {
			L[0][iii0](XDIM) += M[0][iii1]() * r3inv * dX[XDIM];
			L[0][iii0](YDIM) += M[0][iii1]() * r3inv * dX[YDIM];
			L[0][iii0](ZDIM) += M[0][iii1]() * r3inv * dX[ZDIM];
		}

	}

//	printf( "%e\n",sw.get());
}

void grid::compute_ilist() {
	integer lev = nlevel - 2;
	npair np;
	dpair dp;
	for (integer inx = 4; inx <= INX; inx <<= 1) {
		const integer nx = inx + 2 * HBW;
		for (integer i0 = HBW; i0 != nx - HBW; ++i0) {
			for (integer j0 = HBW; j0 != nx - HBW; ++j0) {
				for (integer k0 = HBW; k0 != nx - HBW; ++k0) {
					const integer iii0 = i0 * nx * nx + j0 * nx + k0;
					const integer imin = std::max(integer(HBW), 2 * ((i0 / 2) - 1));
					const integer jmin = std::max(integer(HBW), 2 * ((j0 / 2) - 1));
					const integer kmin = std::max(integer(HBW), 2 * ((k0 / 2) - 1));
					const integer imax = std::min(integer(nx - HBW - 1), 2 * ((i0 / 2) + 1) + 1);
					const integer jmax = std::min(integer(nx - HBW - 1), 2 * ((j0 / 2) + 1) + 1);
					const integer kmax = std::min(integer(nx - HBW - 1), 2 * ((k0 / 2) + 1) + 1);
					for (integer i1 = imin; i1 <= imax; ++i1) {
						for (integer j1 = jmin; j1 <= jmax; ++j1) {
							for (integer k1 = kmin; k1 <= kmax; ++k1) {
								const integer iii1 = i1 * nx * nx + j1 * nx + k1;
								integer max_dist = std::max(std::abs(k0 - k1),
										std::max(std::abs(i0 - i1), std::abs(j0 - j1)));
								if (max_dist > 1 && lev != 0) {
									if (iii1 > iii0) {
										np.lev = lev;
										np.loc.first = iii0;
										np.loc.second = iii1;
										ilist_n.push_back(np);
									}
								} else if (lev == 0 && max_dist > 0) {
									dp.first = iii0;
									dp.second = iii1;
									ilist_d.push_back(dp);
								}
							}
						}
					}
				}
			}
		}
		--lev;
	}
}

void grid::compute_expansions(gsolve_type type) {
	integer lev = nlevel - 1;
	for (integer inx = 4; inx <= INX; inx <<= 1) {
		--lev;
		const integer nxp = (inx / 2) + 2 * HBW;
		const integer nxc = inx + 2 * HBW;

		auto child_index = [=](integer ip, integer jp, integer kp, integer ci) {
			const integer ic = (2 * ip - HBW) + ((ci >> 0) & 1);
			const integer jc = (2 * jp - HBW) + ((ci >> 1) & 1);
			const integer kc = (2 * kp - HBW) + ((ci >> 2) & 1);
			return nxc * nxc * ic + nxc * jc + kc;
		};

		for (integer ip = HBW; ip != nxp - HBW; ++ip) {
			for (integer jp = HBW; jp != nxp - HBW; ++jp) {
				for (integer kp = HBW; kp != nxp - HBW; ++kp) {
					const integer iiip = nxp * nxp * ip + nxp * jp + kp;
					for (integer ci = 0; ci != NVERTEX; ++ci) {
						const integer iiic = child_index(ip, jp, kp, ci);
						const space_vector dX = com[lev][iiic] - com[lev + 1][iiip];
						if (type == RHO) {
							L_c[lev][iiic] += L_c[lev + 1][iiip] << dX;
						}
						L[lev][iiic] += L[lev + 1][iiip] << dX;

					}
				}
			}
		}
	}
}

void grid::compute_multipoles(gsolve_type type) {
	integer lev = 0;
	for (integer inx = INX; inx > 1; inx >>= 1) {

		const integer nxp = inx + 2 * HBW;
		const integer nxc = (2 * inx) + 2 * HBW;

		auto child_index = [=](integer ip, integer jp, integer kp, integer ci) {
			const integer ic = (2 * ip - HBW) + ((ci >> 0) & 1);
			const integer jc = (2 * jp - HBW) + ((ci >> 1) & 1);
			const integer kc = (2 * kp - HBW) + ((ci >> 2) & 1);
			return nxc * nxc * ic + nxc * jc + kc;
		};

		for (integer ip = HBW; ip != nxp - HBW; ++ip) {
			for (integer jp = HBW; jp != nxp - HBW; ++jp) {
				for (integer kp = HBW; kp != nxp - HBW; ++kp) {
					const integer iiip = nxp * nxp * ip + nxp * jp + kp;
					M[lev][iiip] = ZERO;
					if (lev != 0) {
						if (type == RHO) {
							com[lev][iiip] = ZERO;
							real mtot = ZERO;
							for (integer ci = 0; ci != NVERTEX; ++ci) {
								const integer iiic = child_index(ip, jp, kp, ci);
								const real mc = M[lev - 1][iiic]();
								mtot += mc;
								com[lev][iiip] += com[lev - 1][iiic] * mc;
							}
							com[lev][iiip] /= mtot;
						}
						for (integer ci = 0; ci != NVERTEX; ++ci) {
							const integer iiic = child_index(ip, jp, kp, ci);
							const space_vector dX = com[lev - 1][iiic] - com[lev][iiip];
							M[lev][iiip] += M[lev - 1][iiic] >> dX;
						}
					} else {
						if (type == RHO) {
							M[lev][iiip]() = U[rho_i][iiip] * dx * dx * dx;
						} else {
							M[lev][iiip]() = dUdt[rho_i][iiip] * dx * dx * dx;
						}
					}
				}
			}
		}
		++lev;
	}
}
