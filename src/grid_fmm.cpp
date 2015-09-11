/*
 * grid_fmm.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: dmarce1
 */
#include "grid.hpp"
#include "simd.hpp"
#include "stop_watch.hpp"

//#include <boost/thread/lock_guard.hpp>

void grid::solve_gravity(gsolve_type type) {

	compute_multipoles(type);
	compute_interactions(type);
	compute_expansions(type);
}

void grid::compute_interactions(gsolve_type type) {

	npair np;
	dpair dp;
	std::array < simd_vector, NDIM > X;
	std::array<simd_vector, NDIM> Y;
	for (integer lev = 0; lev != nlevel; ++lev) {
		std::fill(L[lev].begin(), L[lev].end(), ZERO);
		std::fill(L_c[lev].begin(), L_c[lev].end(), ZERO);
	}
	const integer list_size = ilist_n.size();
	taylor<4, simd_vector> m0;
	taylor<4, simd_vector> m1;
	taylor<4, simd_vector> n0;
	taylor<4, simd_vector> n1;
	for (integer li = 0; li < list_size; li += simd_len) {
		for (integer i = 0; i != simd_len && li + i < list_size; ++i) {
			const integer iii0 = ilist_n[li + i].loc.first;
			const integer iii1 = ilist_n[li + i].loc.second;
			const integer lev = ilist_n[li + i].lev;
			for (integer d = 0; d != NDIM; ++d) {
				X[d][i] = com[lev][iii0][d];
				Y[d][i] = com[lev][iii1][d];
			}
			for (integer j = 0; j != 20; ++j) {
				m0.ptr()[j][i] = M[lev][iii1].ptr()[j];
				m1.ptr()[j][i] = M[lev][iii0].ptr()[j];
				if (j >= 10) {
					if (type == RHO) {
						n0.ptr()[j][i] = M[lev][iii1].ptr()[j]
								- M[lev][iii0].ptr()[j] * (M[lev][iii1]() / M[lev][iii0]());
					} else {
						n0.ptr()[j][i] = ZERO;
					}
					if (type == RHO) {
						n1.ptr()[j][i] = M[lev][iii0].ptr()[j]
								- M[lev][iii1].ptr()[j] * (M[lev][iii0]() / M[lev][iii1]());
					} else {
						n1.ptr()[j][i] = ZERO;
					}
				}
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
			if (type != RHO) {
				A0() -= m0(a) * D(a);
				A1() += m1(a) * D(a);
			}
			for (integer b = 0; b != NDIM; ++b) {
				const auto tmp = D(a, b) * (real(1) / real(2));
				A0() += m0(a, b) * tmp;
				A1() += m1(a, b) * tmp;
				for (integer c = 0; c != NDIM; ++c) {
					const auto tmp = D(a, b, c) * (real(1) / real(6));
					A0() -= m0(a, b, c) * tmp;
					A1() += m1(a, b, c) * tmp;
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			A0(a) = +m0() * D(a);
			A1(a) = -m1() * D(a);
			for (integer b = 0; b != NDIM; ++b) {
				if (type != RHO) {
					A0(a) -= m0(a) * D(a, b);
					A1(a) -= m1(a) * D(a, b);
				}
				for (integer c = 0; c != NDIM; ++c) {
					const auto tmp = D(a, b, c) * (real(1) / real(2));
					A0(a) += m0(c, b) * tmp;
					A1(a) -= m1(c, b) * tmp;
					if (type == RHO) {
						for (integer d = 0; d != NDIM; ++d) {
							const auto tmp = D(a, b, c, d) * (real(1) / real(6));
							B0(a) -= n0(b, c, d) * tmp;
							B1(a) -= n1(b, c, d) * tmp;
						}
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
		for (integer i = 0; i != simd_len && i + li < list_size; ++i) {
			const integer iii0 = ilist_n[li + i].loc.first;
			const integer iii1 = ilist_n[li + i].loc.second;
			const integer lev = ilist_n[li + i].lev;
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

	const integer dsize = ilist_d.size();
	const integer lev = 0;
	for (integer li = 0; li < dsize; li += simd_len) {
		simd_vector m0, m1;
		for (integer i = 0; i != simd_len && li + i < dsize; ++i) {
			const integer iii0 = ilist_d[li + i].first;
			const integer iii1 = ilist_d[li + i].second;
			for (integer d = 0; d != NDIM; ++d) {
				X[d][i] = com[lev][iii0][d];
				Y[d][i] = com[lev][iii1][d];
			}
			m0[i] = M[lev][iii1]();
			m1[i] = M[lev][iii0]();
		}
		simd_vector phi0, phi1, gx0, gx1, gy0, gy1, gz0, gz1;
		std::array<simd_vector, NDIM> dX;
		simd_vector r = ZERO;
		for (integer d = 0; d != NDIM; ++d) {
			dX[d] = X[d] - Y[d];
			r += dX[d] * dX[d];
		}
		r = sqrt(r);
		const simd_vector rinv = ONE / r;
		const simd_vector r3inv = ONE / (r * r * r);
		phi0 = -m0 * rinv;
		phi1 = -m1 * rinv;
		for (integer d = 0; d != NDIM; ++d) {
			dX[d] *= r3inv;
		}
		if (type == RHO) {
			gx0 = +m0 * dX[XDIM];
			gy0 = +m0 * dX[YDIM];
			gz0 = +m0 * dX[ZDIM];
			gx1 = -m1 * dX[XDIM];
			gy1 = -m1 * dX[YDIM];
			gz1 = -m1 * dX[ZDIM];
		}
		for (integer i = 0; i != simd_len && i + li < dsize; ++i) {
			const integer iii0 = ilist_d[li + i].first;
			const integer iii1 = ilist_d[li + i].second;
			L[lev][iii0]() += phi0[i];
			L[lev][iii1]() += phi1[i];
			if( type == RHO ) {
				L[lev][iii1](XDIM) += gx1[i];
				L[lev][iii1](YDIM) += gy1[i];
				L[lev][iii1](ZDIM) += gz1[i];
				L[lev][iii0](XDIM) += gx0[i];
				L[lev][iii0](YDIM) += gy0[i];
				L[lev][iii0](ZDIM) += gz0[i];
			}
		}
	}
}

void grid::compute_boundary_interactions(gsolve_type type, integer face) {

	npair np;
	dpair dp;
	std::array < simd_vector, NDIM > X;
	std::array<simd_vector, NDIM> Y;
	const integer list_size = ilist_n_bnd[face].size();
	for (integer li = 0; li < list_size; li += simd_len) {
		taylor<4, simd_vector> m0;
		taylor<4, simd_vector> n0;
		for (integer i = 0; i != simd_len && li + i < list_size; ++i) {
			const integer iii0 = ilist_n_bnd[face][li + i].loc.first;
			const integer iii1 = ilist_n_bnd[face][li + i].loc.second;
			const integer lev = ilist_n_bnd[face][li + i].lev;
			for (integer d = 0; d != NDIM; ++d) {
				X[d][i] = com[lev][iii0][d];
				Y[d][i] = com[lev][iii1][d];
			}
			for (integer j = 0; j != 20; ++j) {
				m0.ptr()[j][i] = M[lev][iii1].ptr()[j];
				if (j >= 10) {
					if (type == RHO) {
						n0.ptr()[j][i] = M[lev][iii1].ptr()[j]
								- M[lev][iii0].ptr()[j] * (M[lev][iii1]() / M[lev][iii0]());
					} else {
						n0.ptr()[j][i] = ZERO;
					}
				}
			}
		}
		std::array<simd_vector, NDIM> dX;
		for (integer d = 0; d != NDIM; ++d) {
			dX[d] = X[d] - Y[d];
		}
		taylor<5, simd_vector> D;
		taylor<4, simd_vector> A0, B0;
		D.set_basis(dX);
		B0 = ZERO;
		A0() = m0() * D();
		for (integer a = 0; a != NDIM; ++a) {
			if (type != RHO) {
				A0() -= m0(a) * D(a);
			}
			for (integer b = 0; b != NDIM; ++b) {
				const auto tmp = D(a, b) * (real(1) / real(2));
				A0() += m0(a, b) * tmp;
				for (integer c = 0; c != NDIM; ++c) {
					const auto tmp = D(a, b, c) * (real(1) / real(6));
					A0() -= m0(a, b, c) * tmp;
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			A0(a) = +m0() * D(a);
			for (integer b = 0; b != NDIM; ++b) {
				if (type != RHO) {
					A0(a) -= m0(a) * D(a, b);
				}
				for (integer c = 0; c != NDIM; ++c) {
					const auto tmp = D(a, b, c) * (real(1) / real(2));
					A0(a) += m0(c, b) * tmp;
					if (type == RHO) {
						for (integer d = 0; d != NDIM; ++d) {
							const auto tmp = D(a, b, c, d) * (real(1) / real(6));
							B0(a) -= n0(b, c, d) * tmp;
						}
					}
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			for (integer b = a; b != NDIM; ++b) {
				A0(a, b) = m0() * D(a, b);
				for (integer c = 0; c != NDIM; ++c) {
					A0(a, b) -= m0(c) * D(a, b, c);
				}
			}
		}

		for (integer a = 0; a != NDIM; ++a) {
			for (integer b = a; b != NDIM; ++b) {
				for (integer c = b; c != NDIM; ++c) {
					A0(a, b, c) = +m0() * D(a, b, c);
				}
			}
		}
		for (integer i = 0; i != simd_len && i + li < list_size; ++i) {
			const integer iii0 = ilist_n_bnd[face][li + i].loc.first;
			const integer lev = ilist_n_bnd[face][li + i].lev;
			for (integer j = 0; j != 20; ++j) {
				L[lev][iii0].ptr()[j] += A0.ptr()[j][i];
				if (type == RHO) {
					L_c[lev][iii0].ptr()[j] += B0.ptr()[j][i];
				}
			}
		}
	}
	const integer dsize = ilist_d_bnd[face].size();
	const integer lev = 0;
	for (integer li = 0; li < dsize; li += simd_len) {
		simd_vector m0;
		for (integer i = 0; i != simd_len && li + i < dsize; ++i) {
			const integer iii0 = ilist_d_bnd[face][li + i].first;
			const integer iii1 = ilist_d_bnd[face][li + i].second;
			for (integer d = 0; d != NDIM; ++d) {
				X[d][i] = com[lev][iii0][d];
				Y[d][i] = com[lev][iii1][d];
			}
			m0[i] = M[lev][iii1]();
		}
		simd_vector phi0, gx0, gy0, gz0;
		std::array<simd_vector, NDIM> dX;
		simd_vector r = ZERO;
		for (integer d = 0; d != NDIM; ++d) {
			dX[d] = X[d] - Y[d];
			r += dX[d] * dX[d];
		}
		r = sqrt(r);
		const simd_vector rinv = ONE / r;
		const simd_vector r3inv = ONE / (r * r * r);
		phi0 = -m0 * rinv;
		for (integer d = 0; d != NDIM; ++d) {
			dX[d] *= r3inv;
		}
		if (type == RHO) {
			gx0 = +m0 * dX[XDIM];
			gy0 = +m0 * dX[YDIM];
			gz0 = +m0 * dX[ZDIM];
		}
		for (integer i = 0; i != simd_len && i + li < dsize; ++i) {
			const integer iii0 = ilist_d_bnd[face][li + i].first;
			L[lev][iii0]() += phi0[i];
			if (type == RHO) {
				L[lev][iii0](XDIM) += gx0[i];
				L[lev][iii0](YDIM) += gy0[i];
				L[lev][iii0](ZDIM) += gz0[i];
			}
		}
	}
}

void grid::compute_ilist() {
	std::vector<std::vector<integer> > face_num(nlevel, std::vector<integer>(HN3, integer(0)));
	integer lev = nlevel - 2;
	for (integer inx = 4; inx <= INX; inx <<= 1) {
		const integer nx = inx + 2 * HBW;
		for (integer i0 = 0; i0 != nx; ++i0) {
			for (integer j0 = 0; j0 != nx; ++j0) {
				for (integer k0 = 0; k0 != nx; ++k0) {
					const integer iii0 = i0 * nx * nx + j0 * nx + k0;
					auto& ref = face_num[lev][iii0];
					if (k0 < HBW) {
						ref = FZM;
					} else if (k0 >= nx - HBW) {
						ref = FZP;
					} else if (j0 < HBW) {
						ref = FYM;
					} else if (j0 >= nx - HBW) {
						ref = FYP;
					} else if (i0 < HBW) {
						ref = FXM;
					} else if (i0 >= nx - HBW) {
						ref = FXP;
					} else {
						ref = -1;
					}
				}
			}
		}
		--lev;
	}
	lev = nlevel - 2;
	npair np;
	dpair dp;
	std::vector<npair> ilist_n0;
	std::vector<dpair> ilist_d0;
	std::array<std::vector<npair>, NFACE> ilist_n0_bnd;
	std::array<std::vector<dpair>, NFACE> ilist_d0_bnd;
	for (integer inx = 4; inx <= INX; inx <<= 1) {
		if (is_root || lev == 0) {
			const integer nx = inx + 2 * HBW;
			for (integer i0 = 0; i0 != nx; ++i0) {
				for (integer j0 = 0; j0 != nx; ++j0) {
					for (integer k0 = 0; k0 != nx; ++k0) {
						const integer iii0 = i0 * nx * nx + j0 * nx + k0;
						const integer imin = std::max(integer(0), 2 * ((i0 / 2) - 1));
						const integer jmin = std::max(integer(0), 2 * ((j0 / 2) - 1));
						const integer kmin = std::max(integer(0), 2 * ((k0 / 2) - 1));
						const integer imax = std::min(integer(nx - 1), 2 * ((i0 / 2) + 1) + 1);
						const integer jmax = std::min(integer(nx - 1), 2 * ((j0 / 2) + 1) + 1);
						const integer kmax = std::min(integer(nx - 1), 2 * ((k0 / 2) + 1) + 1);
						for (integer i1 = imin; i1 <= imax; ++i1) {
							for (integer j1 = jmin; j1 <= jmax; ++j1) {
								for (integer k1 = kmin; k1 <= kmax; ++k1) {
									const integer iii1 = i1 * nx * nx + j1 * nx + k1;
									integer max_dist = std::max(std::abs(k0 - k1),
											std::max(std::abs(i0 - i1), std::abs(j0 - j1)));
									if ((max_dist > 1) && ((lev != 0) || !is_leaf)) {
										np.lev = lev;
										np.loc.first = iii0;
										np.loc.second = iii1;
										if (face_num[lev][iii1] == -1 && face_num[lev][iii0] == -1) {
											if (iii1 > iii0) {
												ilist_n0.push_back(np);
											}
										} else if (face_num[lev][iii0] == -1) {
											ilist_n0_bnd[face_num[lev][iii1]].push_back(np);
										}
									} else if ((lev == 0) && (max_dist > 0) && is_leaf) {
										dp.first = iii0;
										dp.second = iii1;
										if (face_num[lev][iii1] == -1 && face_num[lev][iii0] == -1) {
											if (iii1 > iii0) {
												ilist_d0.push_back(dp);
											}
										} else if (face_num[lev][iii0] == -1) {
											ilist_d0_bnd[face_num[lev][iii1]].push_back(dp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		--lev;
	}
	ilist_n = std::vector<npair>(ilist_n0.begin(), ilist_n0.end());
	ilist_d = std::vector<dpair>(ilist_d0.begin(), ilist_d0.end());
	for (integer face = 0; face != NFACE; ++face) {
		ilist_n_bnd[face] = std::vector<npair>(ilist_n0_bnd[face].begin(), ilist_n0_bnd[face].end());
		ilist_d_bnd[face] = std::vector<dpair>(ilist_d0_bnd[face].begin(), ilist_d0_bnd[face].end());
	}
}

expansion_pass_type grid::compute_expansions(gsolve_type type, const expansion_pass_type* parent_expansions) {
	integer lev = is_root ? nlevel - 1 : 1;
	expansion_pass_type exp_ret;
	if (!is_leaf) {
		exp_ret.first.resize(INX * INX * INX);
		if (type == RHO) {
			exp_ret.second.resize(INX * INX * INX);
		}
	}
	for (integer inx = is_root ? 4 : INX; inx <= INX; inx <<= 1) {
		--lev;
		const integer nxp = (inx / 2) + 2 * HBW;

		auto child_index = [=](integer ip, integer jp, integer kp, integer ci, integer bw=HBW) -> integer {
			const integer ic = (2 * (ip - HBW)+bw) + ((ci >> 0) & 1);
			const integer jc = (2 * (jp - HBW)+bw) + ((ci >> 1) & 1);
			const integer kc = (2 * (kp - HBW)+bw) + ((ci >> 2) & 1);
			return (inx+2*bw) * (inx+2*bw) * ic + (inx+2*bw) * jc + kc;
		};

		for (integer ip = HBW; ip != nxp - HBW; ++ip) {
			for (integer jp = HBW; jp != nxp - HBW; ++jp) {
				for (integer kp = HBW; kp != nxp - HBW; ++kp) {
					const integer iiip = nxp * nxp * ip + nxp * jp + kp;
					if (!((lev == 1) && !is_root)) {
						std::array < simd_vector, NDIM > X;
						std::array<simd_vector, NDIM> dX;
						taylor<4, simd_vector> l, lc;
						if (!is_root && lev == 0) {
							const integer index = (INX * INX / 4) * (ip - HBW) + (INX / 2) * (jp - HBW) + (kp - HBW);
							for (integer j = 0; j != 20; ++j) {
								l.ptr()[j] = parent_expansions->first[index].ptr()[j];
								if (type == RHO) {
									lc.ptr()[j] = parent_expansions->second[index].ptr()[j];
								}
							}
						} else {
							for (integer j = 0; j != 20; ++j) {
								l.ptr()[j] = L[lev + 1][iiip].ptr()[j];
								lc.ptr()[j] = L_c[lev + 1][iiip].ptr()[j];
							}
						}
						for (integer ci = 0; ci != NCHILD; ++ci) {
							const integer iiic = child_index(ip, jp, kp, ci);
							for (integer d = 0; d != NDIM; ++d) {
								X[d][ci] = com[lev][iiic][d];
							}
						}
						const auto& Y = com[lev + 1][iiip];
						for (integer d = 0; d != NDIM; ++d) {
							dX[d] = X[d] - Y[d];
						}
						l <<= dX;
						lc <<= dX;
						for (integer ci = 0; ci != NCHILD; ++ci) {
							const integer iiic = child_index(ip, jp, kp, ci);
							for (integer j = 0; j != 20; ++j) {
								if (type == RHO) {
									L_c[lev][iiic].ptr()[j] += lc.ptr()[j][ci];
								}
								L[lev][iiic].ptr()[j] += l.ptr()[j][ci];
							}

							if (!is_leaf && lev == 0) {
								integer index = child_index(ip, jp, kp, ci, 0);
								exp_ret.first[index] = L[lev][iiic];
								if (type == RHO) {
									exp_ret.second[index] = L_c[lev][iiic];
								}
							}
						}
					}
				}
			}
		}
	}

	if (is_leaf) {
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
	return exp_ret;
}

multipole_pass_type grid::compute_multipoles(gsolve_type type, const multipole_pass_type* child_poles) {
	integer lev = 0;
	const real dx3 = dx * dx * dx;
	multipole_pass_type mret;
	if (!is_root) {
		mret.first.resize(INX * INX * INX / NCHILD);
		mret.second.resize(INX * INX * INX / NCHILD);
	}
	integer index = 0;
	for (integer inx = INX; (inx >= std::max(INX / 2, integer(1))) || (inx > 1 && is_root); inx >>= 1) {

		const integer nxp = inx + 2 * HBW;
		const integer nxc = (2 * inx) + 2 * HBW;

		auto child_index = [=](integer ip, integer jp, integer kp, integer ci) -> integer {
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
							simd_vector mc;
							std::array < simd_vector, NDIM > X;
							for (integer ci = 0; ci != NCHILD; ++ci) {
								const integer iiic = child_index(ip, jp, kp, ci);
								mc[ci] = M[lev - 1][iiic]();
								for (integer d = 0; d != NDIM; ++d) {
									X[d][ci] = com[lev - 1][iiic][d];
								}
							}
							real mtot = mc.sum();
							for (integer d = 0; d != NDIM; ++d) {
								com[lev][iiip][d] = (X[d] * mc).sum() / mtot;
							}
						}
						taylor<4, simd_vector> mc, mp;
						std::array<simd_vector, NDIM> x, y, dx;
						for (integer ci = 0; ci != NCHILD; ++ci) {
							const integer iiic = child_index(ip, jp, kp, ci);
							const space_vector& X = com[lev - 1][iiic];
							for (integer j = 0; j != 20; ++j) {
								mc.ptr()[j][ci] = M[lev - 1][iiic].ptr()[j];
								for (integer d = 0; d != NDIM; ++d) {
									x[d][ci] = X[d];
								}
							}
						}
						const space_vector& Y = com[lev][iiip];
						for (integer d = 0; d != NDIM; ++d) {
							simd_vector y = Y[d];
							dx[d] = x[d] - y;
						}
						mp = mc >> dx;
						for (integer j = 0; j != 20; ++j) {
							M[lev][iiip].ptr()[j] = mp.ptr()[j].sum();
						}
					} else {
						if (child_poles == nullptr) {
							if (type == RHO) {
								M[lev][iiip]() = U[rho_i][iiip] * dx3;
							} else {
								M[lev][iiip]() = dUdt[rho_i][iiip] * dx3;
							}
						} else {
							M[lev][iiip] = child_poles->first[index];
							if (type == RHO) {
								com[lev][iiip] = child_poles->second[index];
							}
							++index;
						}
					}
					if (!is_root && (lev == 1)) {
						mret.first[index] = M[lev][iiip];
						mret.second[index] = com[lev][iiip];
						++index;
					}
				}
			}
		}
		++lev;
		index = 0;
	}
	return mret;
}
