/*
 * grid_fmm.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: dmarce1
 */
#include "grid.hpp"

void grid::solve_gravity(gsolve_type type) {
	compute_multipoles(type);
	compute_interactions(type);
	compute_expansions(type);
	if (type == RHO) {
		for (integer iii = 0; iii != HN3; ++iii) {
			G[phi_i][iii] = L[phi_i][0][iii]();
			G[gx_i][iii] = -L[gx_i][0][iii]();
			G[gy_i][iii] = -L[gy_i][0][iii]();
			G[gz_i][iii] = -L[gz_i][0][iii]();
			U[pot_i][iii] = G[phi_i][iii] * U[rho_i][iii];
		}
	} else {
		for (integer iii = 0; iii != HN3; ++iii) {
			dphi_dt[iii] = L[phi_i][0][iii]();
		}
	}
}

void grid::compute_interactions(gsolve_type type) {
	npair np;
	dpair dp;
	for (integer lev = 0; lev != nlevel; ++lev) {
		for (integer d = 0; d != NGF; ++d) {
			std::fill(L[d][lev].begin(), L[d][lev].end(), ZERO);
		}
	}
	for (auto iter = ilist_n.begin(); iter != ilist_n.end(); ++iter) {
		const integer iii0 = iter->loc.first;
		const integer iii1 = iter->loc.second;
		const integer lev = iter->lev;
		space_vector dX = (com[lev][iii0] - com[lev][iii1]);
		taylor<5, real> D;
		D.set_basis(dX);
		taylor<4, real> A, B;
		const multipole& m = M[lev][iii1];
		multipole n;
		if (type == RHO) {
			n = m - M[lev][iii0] * (m() / M[lev][iii0]());
		} else {
			n = ZERO;
		}
		B = ZERO;
		A() = m() * D();
		for (integer a = 0; a != NDIM; ++a) {
			A() -= m(a) * D(a);
			for (integer b = 0; b != NDIM; ++b) {
				A() += m(a, b) * D(a, b) / real(2);
				for (integer c = 0; c != NDIM; ++c) {
					A() -= m(a, b, c) * D(a, b, c) / real(6);
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			A(a) = m() * D(a);
			for (integer b = 0; b != NDIM; ++b) {
				A(a) -= m(a) * D(a, b);
				for (integer c = 0; c != NDIM; ++c) {
					A(a) += m(c, b) * D(a, b, c) / real(2);
					for (integer d = 0; d != NDIM; ++d) {
						B(a) -= n(b, c, d) * D(a, b, c, d) / real(6);
					}
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			for (integer b = a; b != NDIM; ++b) {
				A(a, b) = m() * D(a, b);
				for (integer c = 0; c != NDIM; ++c) {
					A(a, b) -= m(c) * D(a, b, c);
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			for (integer b = a; b != NDIM; ++b) {
				for (integer c = b; c != NDIM; ++c) {
					A(a, b, c) = m() * D(a, b, c);
				}
			}
		}

		L[phi_i][lev][iii0] += A;
		if (type == RHO) {
			A += B;
			for (integer a = 0; a != NDIM; ++a) {
				L[gx_i + a][lev][iii0]() += A(a);
				for (integer b = 0; b != NDIM; ++b) {
					L[gx_i + a][lev][iii0](b) += A(a, b);
					for (integer c = b; c != NDIM; ++c) {
						L[gx_i + a][lev][iii0](b, c) += A(a, b, c);
					}
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
		L[phi_i][0][iii0]() -= M[0][iii1]() * rinv;
		if (type == RHO) {
			L[gx_i][0][iii0]() += M[0][iii1]() * r3inv * dX[XDIM];
			L[gy_i][0][iii0]() += M[0][iii1]() * r3inv * dX[YDIM];
			L[gz_i][0][iii0]() += M[0][iii1]() * r3inv * dX[ZDIM];
		}

	}
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
									np.lev = lev;
									np.loc.first = iii0;
									np.loc.second = iii1;
									ilist_n.push_back(np);
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
							for (integer d = 0; d != NGF; ++d) {
								L[d][lev][iiic] += L[d][lev + 1][iiip] << dX;
							}
						} else {
							L[phi_i][lev][iiic] += L[phi_i][lev + 1][iiip] << dX;
						}
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
