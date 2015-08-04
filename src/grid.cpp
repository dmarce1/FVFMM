/*
 * grid.cpp
 *
 *  Created on: May 26, 2015
 *      Author: dmarce1
 */

#include "grid.hpp"
#include <cmath>
#include <cassert>

real grid::omega = DEFAULT_OMEGA;
//real grid::omega = ZERO;

std::vector<real> grid::conserved_sums() const {
	std::vector<real> sum(NF);
	const real dV = dx*dx*dx;
	std::fill(sum.begin(), sum.end(), ZERO);
	for( integer i = HBW; i != HNX - HBW; ++i) {
		for( integer j= HBW; j != HNX -HBW; ++j) {
			for( integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				for( integer field = 0; field != NF; ++field) {
					sum[field] += U[field][iii] * dV;
					if( field == egas_i ) {
						sum[field] += U[pot_i][iii] * HALF * dV;
					}
				}
			}
		}
	}
	return sum;
}

std::vector<real> grid::s_sums() const {
	std::vector<real> sum(NDIM);
	const real dV = dx*dx*dx;
	std::fill(sum.begin(), sum.end(), ZERO);
	for( integer i = HBW; i != HNX - HBW; ++i) {
		for( integer j= HBW; j != HNX -HBW; ++j) {
			for( integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				for( integer d = 0; d != NDIM; ++d ) {
					sum[d] += S[d][iii] * dV;
//					sum[d] += omega * dx * dx * U[rho_i][iii];
				}
			}
		}
	}
	for( integer d = 0; d != NDIM; ++d) {
		sum[d] -= S_out[d];
	}
	return sum;
}

std::vector<real> grid::l_sums() const {
	std::vector<real> sum(NDIM);
	const real dV = dx*dx*dx;
	std::fill(sum.begin(), sum.end(), ZERO);
	for( integer i = HBW; i != HNX - HBW; ++i) {
		for( integer j= HBW; j != HNX -HBW; ++j) {
			for( integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				sum[XDIM] += X[YDIM][iii] * U[sz_i][iii] * dV;
				sum[XDIM] -= X[ZDIM][iii] * U[sy_i][iii] * dV;

				sum[YDIM] -= X[XDIM][iii] * U[sz_i][iii] * dV;
				sum[YDIM] += X[ZDIM][iii] * U[sx_i][iii] * dV;

				sum[ZDIM] += X[XDIM][iii] * U[sy_i][iii] * dV;
				sum[ZDIM] -= X[YDIM][iii] * U[sx_i][iii] * dV;

			}
		}
	}
	return sum;
}

std::vector<real> grid::conserved_outflows() const {
	auto Uret = U_out;
	Uret[egas_i] += Uret[pot_i];
	return Uret;
}

void grid::set_omega(real o) {
	omega = o;
}

real grid::get_omega() {
	return omega;
}


bool grid::refine_me(integer lev) const {
	return lev < MAX_LEVEL;
}

integer grid::level_count() const {
	return nlevel;
}

grid::~grid() {

}


real& grid::spin_value(integer f, integer i, integer j, integer k) {
	return S[f][i * DNX + j * DNY + k * DNZ];
}


real grid::spin_value(integer f, integer i, integer j, integer k) const {
	return S[f][i * DNX + j * DNY + k * DNZ];
}


real& grid::hydro_value(integer f, integer i, integer j, integer k) {
	return U[f][i * DNX + j * DNY + k * DNZ];
}


space_vector& grid::center_of_mass_value(integer i, integer j, integer k) {
	return com[0][i * HNX * HNX + j * HNX + k];
}


const space_vector& grid::center_of_mass_value(integer i, integer j, integer k) const {
	return com[0][i * HNX * HNX + j * HNX + k];
}


multipole& grid::multipole_value(integer lev, integer i, integer j, integer k) {
	const integer bw = HBW;
	const integer inx = INX >> lev;
	const integer nx = 2 * bw + inx;
	return M[lev][i * nx * nx + j * nx + k];
}

const multipole& grid::multipole_value(integer lev, integer i, integer j, integer k) const {
	const integer bw = HBW;
	const integer inx = INX >> lev;
	const integer nx = 2 * bw + inx;
	return M[lev][i * nx * nx + j * nx + k];
}

real grid::hydro_value(integer f, integer i, integer j, integer k) const {
	return U[f][i * DNX + j * DNY + k * DNZ];
}

inline real minmod(real a, real b) {
	return (std::copysign(HALF, a) + std::copysign(HALF, b)) * std::min(std::abs(a), std::abs(b));
}

inline real minmod_theta(real a, real b, real theta) {
	return minmod(theta * minmod(a, b), HALF * (a + b));
}

grid::grid(real _dx, std::array<real, NDIM> _xmin, integer flags) : U(NF), U0(NF), dUdt(NF), Uf(NFACE), F(NDIM), X(NDIM), G(NGF), G0(NGF), S0(NDIM), S(NDIM), src(NF), ilist_d_bnd(NFACE),ilist_n_bnd(NFACE),
		is_root(flags & GRID_IS_ROOT), is_leaf(flags & GRID_IS_LEAF) {
	dx = _dx;
	xmin = _xmin;
	allocate();
}

void grid::allocate() {
	t = ZERO;
	step_num = 0;
	U_out0 = std::vector<real>(NF, ZERO);
	U_out = std::vector<real>(NF, ZERO);
	S_out0 = std::vector<real>(NF, ZERO);
	S_out = std::vector<real>(NF, ZERO);
	dphi_dt = std::vector<real>(HN3);
	for (integer field = 0; field != NGF; ++field) {
		G[field].resize(HN3);
		G0[field].resize(HN3);
	}
	for (integer dim = 0; dim != NDIM; ++dim) {
		X[dim].resize(HN3);
	}
	for (integer i = 0; i != HNX; ++i) {
		for (integer j = 0; j != HNX; ++j) {
			for (integer k = 0; k != HNX; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				X[XDIM][iii] = (real(i - HBW) + HALF) * dx + xmin[XDIM];
				X[YDIM][iii] = (real(j - HBW) + HALF) * dx + xmin[YDIM];
				X[ZDIM][iii] = (real(k - HBW) + HALF) * dx + xmin[ZDIM];
			}
		}
	}
	for (integer field = 0; field != NDIM; ++field) {
		S0[field].resize(HN3);
		S[field].resize(HN3);
	}
	for (integer field = 0; field != NF; ++field) {
		src[field].resize(HN3);
		U0[field].resize(HN3);
		U[field].resize(HN3);
		dUdt[field].resize(HN3);
		for (integer dim = 0; dim != NDIM; ++dim) {
			F[dim][field].resize(HN3);
		}
		for (integer face = 0; face != NFACE; ++face) {
			Uf[face][field].resize(HN3);
		}
	}
	nlevel = 0;
	for (integer inx = INX; inx > 1; inx /= 2) {
		++nlevel;
	}
	com.resize(nlevel);
	M.resize(nlevel);
	L.resize(nlevel);
	L_c.resize(nlevel);
	nlevel = 0;
	for (integer inx = INX; inx > 1; inx /= 2) {
		const integer this_nx = inx + 2 * HBW;
		const integer sz = this_nx * this_nx * this_nx;
		com[nlevel].resize(sz);
		if (is_root || (inx >= INX / 2)) {
			M[nlevel].resize(sz);
		}
		if (is_root || (inx >= INX)) {
			L[nlevel].resize(sz);
			L_c[nlevel].resize(sz);
		}
		++nlevel;
	}
	for (integer iii = 0; iii != HN3; ++iii) {
		for (integer dim = 0; dim != NDIM; ++dim) {
			com[0][iii][dim] = X[dim][iii];
		}
	}
	compute_ilist();

}

grid::grid() : U(NF), U0(NF), dUdt(NF), Uf(NFACE), F(NDIM), X(NDIM), G(NGF), G0(NGF), S0(NDIM), S(NDIM), src(NF), ilist_d_bnd(NFACE), ilist_n_bnd(NFACE)
{
}

grid::grid(const std::function<std::vector<real>(real, real, real)>& init_func, real _dx, std::array<real, NDIM> _xmin,
		integer flags) : U(NF), U0(NF), dUdt(NF), Uf(NFACE), F(NDIM), X(NDIM), G(NGF), G0(NGF), S0(NDIM), S(NDIM), src(NF), ilist_d_bnd(NFACE),ilist_n_bnd(NFACE),
		is_root(flags & GRID_IS_ROOT), is_leaf(flags & GRID_IS_LEAF), U_out(NF, ZERO), U_out0(NF, ZERO), S_out(NDIM,
				ZERO), S_out0(NDIM, ZERO), dphi_dt(HN3) {
	dx = _dx;
	xmin = _xmin;
	allocate();
	for (integer i = 0; i != HNX; ++i) {
		for (integer j = 0; j != HNX; ++j) {
			for (integer k = 0; k != HNX; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				std::vector<real> this_u = init_func(X[XDIM][iii], X[YDIM][iii], X[ZDIM][iii]);
				for (integer field = 0; field != NF; ++field) {
					U[field][iii] = this_u[field];
				}
				for (integer d = 0; d != NDIM; ++d) {
					S[d][iii] = ZERO;

				}
			}
		}
	}
}

void grid::reconstruct() {

	std::array<std::vector<real>, NF> slpx, slpy, slpz;
	for (integer field = 0; field != NF; ++field) {
		slpx[field].resize(HN3);
		slpy[field].resize(HN3);
		slpz[field].resize(HN3);
		if (field != rho_i) {
#pragma GCC ivdep
			for (integer iii = 0; iii != HN3; ++iii) {
				U[field][iii] /= U[rho_i][iii];
			}
		}
		if (field == sx_i) {
#pragma GCC ivdep
			for (integer iii = 0; iii != HN3; ++iii) {
				U[field][iii] += omega * X[YDIM][iii];
			}
		} else if (field == sy_i) {
#pragma GCC ivdep
			for (integer iii = 0; iii != HN3; ++iii) {
				U[field][iii] -= omega * X[XDIM][iii];
			}
		}
#pragma GCC ivdep
		for (integer iii = HNX * HNX; iii != HN3 - HNX * HNX; ++iii) {
			const real u0 = U[field][iii];
			slpx[field][iii] = minmod_theta(U[field][iii + DNX] - u0, u0 - U[field][iii - DNX], 1.0);
			slpy[field][iii] = minmod_theta(U[field][iii + DNY] - u0, u0 - U[field][iii - DNY], 1.0);
			slpz[field][iii] = minmod_theta(U[field][iii + DNZ] - u0, u0 - U[field][iii - DNZ], 1.0);
		}
	}
#pragma GCC ivdep
	for (integer iii = HNX * HNX; iii != HN3 - HNX * HNX; ++iii) {
		const real lx_sum = slpy[sz_i][iii] + slpz[sy_i][iii];
		const real ly_sum = slpx[sz_i][iii] + slpz[sx_i][iii];
		const real lz_sum = slpx[sy_i][iii] + slpy[sx_i][iii];
		const real lx_dif = +S[XDIM][iii] / U[rho_i][iii] / dx;
		const real ly_dif = -S[YDIM][iii] / U[rho_i][iii] / dx;
		const real lz_dif = +S[ZDIM][iii] / U[rho_i][iii] / dx/* - omega * dx*/;
		slpx[sy_i][iii] = (lz_sum + lz_dif) * HALF;
		slpy[sx_i][iii] = (lz_sum - lz_dif) * HALF;
		slpx[sz_i][iii] = (ly_sum + ly_dif) * HALF;
		slpz[sx_i][iii] = (ly_sum - ly_dif) * HALF;
		slpy[sz_i][iii] = (lx_sum + lx_dif) * HALF;
		slpz[sy_i][iii] = (lx_sum - lx_dif) * HALF;
	}

	for (integer field = 0; field != NF; ++field) {
#pragma GCC ivdep
		for (integer iii = HNX * HNX; iii != HN3 - HNX * HNX; ++iii) {
			const real u0 = U[field][iii];
			Uf[FXP][field][iii] = u0 + HALF * slpx[field][iii];
			Uf[FXM][field][iii] = u0 - HALF * slpx[field][iii];
			Uf[FYP][field][iii] = u0 + HALF * slpy[field][iii];
			Uf[FYM][field][iii] = u0 - HALF * slpy[field][iii];
			Uf[FZP][field][iii] = u0 + HALF * slpz[field][iii];
			Uf[FZM][field][iii] = u0 - HALF * slpz[field][iii];
		}
		if (field == pot_i) {
#pragma GCC ivdep
			for (integer iii = HNX * HNX; iii != HN3 - HNX * HNX; ++iii) {
				const real phi_x = HALF * (Uf[FXM][field][iii] + Uf[FXP][field][iii - DNX]);
				const real phi_y = HALF * (Uf[FYM][field][iii] + Uf[FYP][field][iii - DNY]);
				const real phi_z = HALF * (Uf[FZM][field][iii] + Uf[FZP][field][iii - DNZ]);
				Uf[FXM][field][iii] = phi_x;
				Uf[FYM][field][iii] = phi_y;
				Uf[FZM][field][iii] = phi_z;
				Uf[FXP][field][iii - DNX] = phi_x;
				Uf[FYP][field][iii - DNY] = phi_y;
				Uf[FZP][field][iii - DNZ] = phi_z;
			}
		}
		else if (field == sx_i) {
#pragma GCC ivdep
			for (integer iii = 0; iii != HN3; ++iii) {
				U[field][iii]       -= omega *  X[YDIM][iii];
			}
#pragma GCC ivdep
			for (integer iii = DNY; iii != HN3-DNY; ++iii) {
				Uf[FXM][field][iii] -= omega *  X[YDIM][iii];
				Uf[FYM][field][iii] -= omega * (X[YDIM][iii] + X[YDIM][iii-DNY])*HALF;
				Uf[FZM][field][iii] -= omega *  X[YDIM][iii];
				Uf[FXP][field][iii] -= omega *  X[YDIM][iii];
				Uf[FYP][field][iii] -= omega * (X[YDIM][iii] + X[YDIM][iii+DNY])*HALF;
				Uf[FZP][field][iii] -= omega *  X[YDIM][iii];
			}
		} else if (field == sy_i) {
#pragma GCC ivdep
			for (integer iii = 0; iii != HN3; ++iii) {
				U[field][iii]       += omega *  X[XDIM][iii];
			}
#pragma GCC ivdep
			for (integer iii = DNX; iii != HN3-DNX; ++iii) {
				Uf[FXM][field][iii] += omega * (X[XDIM][iii] + X[XDIM][iii-DNX])*HALF;
				Uf[FYM][field][iii] += omega *  X[XDIM][iii];
				Uf[FZM][field][iii] += omega *  X[XDIM][iii];
				Uf[FXP][field][iii] += omega * (X[XDIM][iii] + X[XDIM][iii+DNX])*HALF;
				Uf[FYP][field][iii] += omega *  X[XDIM][iii];
				Uf[FZP][field][iii] += omega *  X[XDIM][iii];
			}
		}
		if (field != rho_i) {
#pragma GCC ivdep
			for (integer iii = 0; iii != HN3; ++iii) {
				U[field][iii] *= U[rho_i][iii];
				for (integer face = 0; face != NFACE; ++face) {
					Uf[face][field][iii] *= Uf[face][rho_i][iii];
				}
			}
		}
	}
}

real grid::get_time() const {
	return t;
}

integer grid::get_step() const {
	return step_num;
}

real grid::compute_fluxes() {
	real max_lambda = ZERO;
	std::array<std::vector<real>, NF> ur, ul, f;
	std::vector<space_vector> x;

	const integer line_sz = HNX - 2 * HBW + 1;
	for (integer field = 0; field != NF; ++field) {
		ur[field].resize(line_sz);
		ul[field].resize(line_sz);
		f[field].resize(line_sz);
	}
	x.resize(line_sz);

	for (integer dim = 0; dim != NDIM; ++dim) {

		const integer dx_i = dim;
		const integer dy_i = (dim == XDIM ? YDIM : XDIM);
		const integer dz_i = (dim == ZDIM ? YDIM : ZDIM);
		const integer face_p = 2 * dim + 1;
		const integer face_m = 2 * dim;

		for (integer k = HBW; k != HNX - HBW; ++k) {
			for (integer j = HBW; j != HNX - HBW; ++j) {
				for (integer i = HBW; i != HNX - HBW + 1; ++i) {
					const integer i0 = DN[dx_i] * i + DN[dy_i] * j + DN[dz_i] * k;
					const integer im = i0 - DN[dx_i];
					for (integer field = 0; field != NF; ++field) {
						ur[field][i - HBW] = Uf[face_m][field][i0];
						ul[field][i - HBW] = Uf[face_p][field][im];
					}
					for( integer d = 0; d != NDIM; ++d) {
						x[i - HBW][d] = (X[d][i0] + X[d][im])*HALF;
					}
				}
				const real this_max_lambda = roe_fluxes(f, ul, ur, x, omega, dim);
				max_lambda = std::max(max_lambda, this_max_lambda);
				for (integer field = 0; field != NF; ++field) {
					for (integer i = HBW; i != HNX - HBW + 1; ++i) {
						const integer i0 = DN[dx_i] * i + DN[dy_i] * j + DN[dz_i] * k;
						F[dim][field][i0] = f[field][i - HBW];
			//			if( i == HBW || i == HNX - HBW) {
			//				F[dim][field][i0]  = ZERO;
			//			}
					}
				}
			}
		}
	}
//	if( max_lambda > 1 )
//	printf( "%e\n", max_lambda);

	return max_lambda;
}

void grid::store() {
	for (integer field = 0; field != NF; ++field) {
#pragma GCC ivdep
		for (integer iii = 0; iii != HN3; ++iii) {
			U0[field][iii] = U[field][iii];
		}
	}
	for (integer field = 0; field != NDIM; ++field) {
#pragma GCC ivdep
		for (integer iii = 0; iii != HN3; ++iii) {
			S0[field][iii] = S[field][iii];
		}
	}
	for (integer field = 0; field != NGF; ++field) {
#pragma GCC ivdep
		for (integer iii = 0; iii != HN3; ++iii) {
			G0[field][iii] = G[field][iii];
		}
	}
	U_out0 = U_out;
	S_out0 = S_out;
}

void grid::restore() {
	for (integer field = 0; field != NF; ++field) {
#pragma GCC ivdep
		for (integer iii = 0; iii != HN3; ++iii) {
			U[field][iii] = U0[field][iii];
		}
	}
	for (integer field = 0; field != NDIM; ++field) {
#pragma GCC ivdep
		for (integer iii = 0; iii != HN3; ++iii) {
			S[field][iii] = S0[field][iii];
		}
	}
	for (integer field = 0; field != NGF; ++field) {
#pragma GCC ivdep
		for (integer iii = 0; iii != HN3; ++iii) {
			G[field][iii] = G0[field][iii];
		}
	}
	U_out = U_out0;
	S_out = S_out0;
}

void grid::boundaries() {
	for (integer face = 0; face != NFACE; ++face) {
		set_physical_boundaries(face);
	}
}

void grid::set_physical_boundaries(integer face) {
	const integer dni = face / 2 == XDIM ? DNY : DNX;
	const integer dnj = face / 2 == ZDIM ? DNY : DNZ;
	const integer dnk = DN[face / 2];
	const integer klb = face % 2 == 0 ? 0 : HNX - HBW;
	const integer kub = face % 2 == 0 ? HBW : HNX;
	const integer klb2 = face % 2 == 0 ? 0 : HNX - 2*HBW;
	const integer kub2 = face % 2 == 0 ? 2*HBW : HNX;
	for (integer k = klb2; k != kub2; ++k) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			for (integer i = HBW; i != HNX - HBW; ++i) {
				const integer iii = i * dni + j * dnj + k * dnk;
				U[sx_i][iii] += omega * X[YDIM][iii] * U[rho_i][iii];
				U[sy_i][iii] -= omega * X[XDIM][iii] * U[rho_i][iii];
			}
		}
	}
	for (integer field = 0; field != NF; ++field) {
		for (integer k = klb; k != kub; ++k) {
			for (integer j = HBW; j != HNX - HBW; ++j) {
				for (integer i = HBW; i != HNX - HBW; ++i) {
					integer k0;
					switch( boundary_types[face] == REFLECT) {
					case REFLECT:
						 k0	= face % 2 == 0 ? (2 * HBW - k - 1) : (2 * (HNX - HBW) - k - 1);
						break;
					case OUTFLOW:
						 k0	= face % 2 == 0 ? HBW : HNX - HBW - 1;
						break;
					}
					const real value = U[field][i * dni + j * dnj + k0 * dnk];
					real& ref = U[field][i * dni + j * dnj + k * dnk];
					if (field != sx_i + face / 2) {
						ref = +value;
					} else {
						switch (boundary_types[face]) {
						case REFLECT:
							ref = -value;
							break;
						case OUTFLOW:
							if (face % 2 == 0) {
								ref = std::min(value, ZERO);
							} else {
								ref = std::max(value, ZERO);
							}
							break;
						}
					}
				}
			}
		}
	}
	for (integer k = klb2; k != kub2; ++k) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			for (integer i = HBW; i != HNX - HBW; ++i) {
				const integer iii = i * dni + j * dnj + k * dnk;
				U[sx_i][iii] -= omega * X[YDIM][iii] * U[rho_i][iii];
				U[sy_i][iii] += omega * X[XDIM][iii] * U[rho_i][iii];
			}
		}
	}
}

void grid::compute_sources() {
	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
#pragma GCC ivdep
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				for (integer field = 0; field != NF; ++field) {
					src[field][iii] = ZERO;
				}
				src[sx_i][iii] += U[rho_i][iii] * G[gx_i][iii];
				src[sy_i][iii] += U[rho_i][iii] * G[gy_i][iii];
				src[sz_i][iii] += U[rho_i][iii] * G[gz_i][iii];
				src[sx_i][iii] += omega * U[sy_i][iii];
				src[sy_i][iii] -= omega * U[sx_i][iii];



				src[egas_i][iii] -= omega * X[YDIM][iii] * U[rho_i][iii] * G[gx_i][iii];
				src[egas_i][iii] += omega * X[XDIM][iii] * U[rho_i][iii] * G[gy_i][iii];
			}
		}
	}
}

void grid::compute_dudt() {
	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			for (integer field = 0; field != NF; ++field) {
#pragma GCC ivdep
				for (integer k = HBW; k != HNX - HBW; ++k) {
					const integer iii = DNX * i + DNY * j + DNZ * k;
					dUdt[field][iii] = ZERO;
					dUdt[field][iii] -= (F[XDIM][field][iii + DNX] - F[XDIM][field][iii]) / dx;
					dUdt[field][iii] -= (F[YDIM][field][iii + DNY] - F[YDIM][field][iii]) / dx;
					dUdt[field][iii] -= (F[ZDIM][field][iii + DNZ] - F[ZDIM][field][iii]) / dx;
					dUdt[field][iii] += src[field][iii];
				}
			}
#pragma GCC ivdep
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				dUdt[egas_i][iii] += dUdt[pot_i][iii];
				dUdt[pot_i][iii] = ZERO;
			}
#pragma GCC ivdep
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				dUdt[egas_i][iii] -= (dUdt[rho_i][iii] * G[phi_i][iii]) * HALF;
			}
		}
	}
//	solve_gravity(DRHODT);
}

void grid::egas_to_etot() {
#pragma GCC ivdep
	for (integer iii = 0; iii != HN3; ++iii) {
		U[egas_i][iii] += U[pot_i][iii] * HALF;
	}
}

void grid::etot_to_egas() {
#pragma GCC ivdep
	for (integer iii = 0; iii != HN3; ++iii) {
		U[egas_i][iii] -= U[pot_i][iii] * HALF;
	}
}

void grid::next_u(integer rk, real dt) {

	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
#pragma GCC ivdep
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				dUdt[egas_i][iii] += (dphi_dt[iii] * U[rho_i][iii]) * HALF;
			}
		}
	}

	std::vector<real> ds(NDIM, ZERO);
	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				real d = (-(F[YDIM][sz_i][iii + DNY] + F[YDIM][sz_i][iii])
						+ (F[ZDIM][sy_i][iii + DNZ] + F[ZDIM][sy_i][iii])) * HALF;
				d -= omega * X[ZDIM][iii] * U[sx_i][iii];
				const real s1x = S[XDIM][iii] + d * dt;
				S[XDIM][iii] = (ONE - rk_beta[rk]) * S0[XDIM][iii] + rk_beta[rk] * s1x;
			}
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				real d = (+(F[XDIM][sz_i][iii + DNX] + F[XDIM][sz_i][iii])
						- (F[ZDIM][sx_i][iii + DNZ] + F[ZDIM][sx_i][iii])) * HALF;
				d -= omega * X[ZDIM][iii] * U[sy_i][iii];
				const real s1y = S[YDIM][iii] + d * dt;
				S[YDIM][iii] = (ONE - rk_beta[rk]) * S0[YDIM][iii] + rk_beta[rk] * s1y;
			}
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				real d = (-(F[XDIM][sy_i][iii + DNX] + F[XDIM][sy_i][iii])
						 + (F[YDIM][sx_i][iii + DNY] + F[YDIM][sx_i][iii])) * HALF;
				d += omega * (X[XDIM][iii] * U[sx_i][iii] + X[YDIM][iii] * U[sy_i][iii]);
				const real s1z = S[ZDIM][iii] + d * dt;
				S[ZDIM][iii] = (ONE - rk_beta[rk]) * S0[ZDIM][iii] + rk_beta[rk] * s1z;
			}
			for (integer field = 0; field != NF; ++field) {
				for (integer k = HBW; k != HNX - HBW; ++k) {
					const integer iii = DNX * i + DNY * j + DNZ * k;
					const real u1 = U[field][iii] + dUdt[field][iii] * dt;
					const real u0 = U0[field][iii];
					U[field][iii] = (ONE - rk_beta[rk]) * u0 + rk_beta[rk] * u1;
				}
			}
		}
	}

	std::vector<real> du_out(NF, ZERO);
	std::vector<real> ds_out(NF, ZERO);

	du_out[sx_i] += omega * U_out[sy_i] * dt;
	du_out[sy_i] -= omega * U_out[sx_i] * dt;


	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			const integer iii_p = DNX * (HNX - HBW) + DNY * i + DNZ * j;
			const integer jjj_p = DNY * (HNX - HBW) + DNZ * i + DNX * j;
			const integer kkk_p = DNZ * (HNX - HBW) + DNX * i + DNY * j;
			const integer iii_m = DNX * (HBW) + DNY * i + DNZ * j;
			const integer jjj_m = DNY * (HBW) + DNZ * i + DNX * j;
			const integer kkk_m = DNZ * (HBW) + DNX * i + DNY * j;
			std::vector<real> du(NF);
			for (integer field = 0; field != NF; ++field) {
				du[field] = ZERO;
				du[field] += (F[XDIM][field][iii_p] - F[XDIM][field][iii_m]) * dx * dx;
				du[field] += (F[YDIM][field][jjj_p] - F[YDIM][field][jjj_m]) * dx * dx;
				du[field] += (F[ZDIM][field][kkk_p] - F[ZDIM][field][kkk_m]) * dx * dx;
			}
			//	du[egas_i] += du[pot_i];
			for (integer field = 0; field != NF; ++field) {
				du_out[field] += du[field] * dt;
			}

			const real xp = X[XDIM][iii_p] - HALF * dx;
			const real xm = X[XDIM][iii_m] - HALF * dx;
			const real yp = X[YDIM][jjj_p] - HALF * dx;
			const real ym = X[YDIM][jjj_m] - HALF * dx;
			const real zp = X[ZDIM][kkk_p] - HALF * dx;
			const real zm = X[ZDIM][kkk_m] - HALF * dx;

			ds[XDIM] = -(yp * F[YDIM][sz_i][jjj_p] - ym * F[YDIM][sz_i][jjj_m]) * dx * dx;
			ds[XDIM] += (zp * F[ZDIM][sy_i][kkk_p] - zm * F[ZDIM][sy_i][kkk_m]) * dx * dx;
			ds[XDIM] -= (X[YDIM][iii_p] * F[XDIM][sz_i][iii_p] - X[YDIM][iii_m] * F[XDIM][sz_i][iii_m]) * dx * dx;
			ds[XDIM] -= (X[YDIM][kkk_p] * F[ZDIM][sz_i][kkk_p] - X[YDIM][kkk_m] * F[ZDIM][sz_i][kkk_m]) * dx * dx;
			ds[XDIM] += (X[ZDIM][iii_p] * F[XDIM][sy_i][iii_p] - X[ZDIM][iii_m] * F[XDIM][sy_i][iii_m]) * dx * dx;
			ds[XDIM] += (X[ZDIM][jjj_p] * F[YDIM][sy_i][jjj_p] - X[ZDIM][jjj_m] * F[YDIM][sy_i][jjj_m]) * dx * dx;

			ds[YDIM] = +(xp * F[XDIM][sz_i][iii_p] - xm * F[XDIM][sz_i][iii_m]) * dx * dx;
			ds[YDIM] -= (zp * F[ZDIM][sx_i][kkk_p] - zm * F[ZDIM][sx_i][kkk_m]) * dx * dx;
			ds[YDIM] += (X[XDIM][jjj_p] * F[YDIM][sz_i][jjj_p] - X[XDIM][jjj_m] * F[YDIM][sz_i][jjj_m]) * dx * dx;
			ds[YDIM] += (X[XDIM][kkk_p] * F[ZDIM][sz_i][kkk_p] - X[XDIM][kkk_m] * F[ZDIM][sz_i][kkk_m]) * dx * dx;
			ds[YDIM] -= (X[ZDIM][iii_p] * F[XDIM][sx_i][iii_p] - X[ZDIM][iii_m] * F[XDIM][sx_i][iii_m]) * dx * dx;
			ds[YDIM] -= (X[ZDIM][jjj_p] * F[YDIM][sx_i][jjj_p] - X[ZDIM][jjj_m] * F[YDIM][sx_i][jjj_m]) * dx * dx;

			ds[ZDIM] = -(xp * F[XDIM][sy_i][iii_p] - xm * F[XDIM][sy_i][iii_m]) * dx * dx;
			ds[ZDIM] += (X[YDIM][iii_p] * F[XDIM][sx_i][iii_p] - X[YDIM][iii_m] * F[XDIM][sx_i][iii_m]) * dx * dx;

			ds[ZDIM] -= (X[XDIM][jjj_p] * F[YDIM][sy_i][jjj_p] - X[XDIM][jjj_m] * F[YDIM][sy_i][jjj_m]) * dx * dx;
			ds[ZDIM] += (yp * F[YDIM][sx_i][jjj_p] - ym * F[YDIM][sx_i][jjj_m]) * dx * dx;

			ds[ZDIM] -= (X[XDIM][kkk_p] * F[ZDIM][sy_i][kkk_p] - X[XDIM][kkk_m] * F[ZDIM][sy_i][kkk_m]) * dx * dx;
			ds[ZDIM] += (X[YDIM][kkk_p] * F[ZDIM][sx_i][kkk_p] - X[YDIM][kkk_m] * F[ZDIM][sx_i][kkk_m]) * dx * dx;

			for (integer field = 0; field != NDIM; ++field) {
				ds_out[field] += ds[field] * dt;
			}
		}
	}
#pragma GCC ivdep
	for (integer field = 0; field != NF; ++field) {
		const real out1 = U_out[field] + du_out[field];
		const real out0 = U_out0[field];
		U_out[field] = (ONE - rk_beta[rk]) * out0 + rk_beta[rk] * out1;
	}
	for (integer field = 0; field != NDIM; ++field) {
		const real out1 = S_out[field] + ds_out[field];
		const real out0 = S_out0[field];
		S_out[field] = (ONE - rk_beta[rk]) * out0 + rk_beta[rk] * out1;
	}
	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = DNX * i + DNY * j + DNZ * k;
				real ek = ZERO;
				ek += HALF * std::pow(U[sx_i][iii], 2) / U[rho_i][iii];
				ek += HALF * std::pow(U[sy_i][iii], 2) / U[rho_i][iii];
				ek += HALF * std::pow(U[sz_i][iii], 2) / U[rho_i][iii];
				real ei = U[egas_i][iii] - ek;
				if (ei > de_switch1 * U[egas_i][iii]) {
					U[tau_i][iii] = std::pow(ei, ONE / fgamma);
				}
				if( U[tau_i][iii] < ZERO ) {
					printf( "Tau is negative\n");
					abort();
				} else if( U[rho_i][iii] <= ZERO ) {
					printf( "Rho is non-positive\n");
					abort();
				}
			}
		}
	}
}

real grid::step() {

	bool cfl_redo = false;
	real cfl0 = cfl;
	real dt, a;
	store();
	do {

		reconstruct();
		a = compute_fluxes();
		compute_sources();
		dt = cfl0 * dx / a;
		compute_dudt();
		next_u(0, dt);
		egas_to_etot();
		solve_gravity();
		etot_to_egas();
		boundaries();

		reconstruct();
		a = compute_fluxes();
		compute_sources();
		if (cfl0 * dx / a < HALF * dt) {
			printf("REDO SUBSTEP\n");
			cfl0 *= HALF;
			restore();
			cfl_redo = true;
		} else {
			compute_dudt();
			next_u(1, dt);
			egas_to_etot();
			solve_gravity();
			etot_to_egas();
			boundaries();
			cfl_redo = false;
		}

	} while (cfl_redo);
	t += dt;
	diagnostics();
	++step_num;
	return dt;
}

void grid::diagnostics() {
	std::vector<real> sum(NF);
	real gx_sum = ZERO, gy_sum = ZERO, gz_sum = ZERO;
	real glx_sum = ZERO, gly_sum = ZERO, glz_sum = ZERO;
	real sx_sum = ZERO, sy_sum = ZERO, sz_sum = ZERO;
	real lx_sum = ZERO, ly_sum = ZERO, lz_sum = ZERO;
	real dV = std::pow(dx, NDIM);
	for (integer i = HBW; i != HNX - HBW; ++i) {
		for (integer j = HBW; j != HNX - HBW; ++j) {
			for (integer field = 0; field != NF; ++field) {
				for (integer k = HBW; k != HNX - HBW; ++k) {
					const integer iii = i * DNX + j * DNY + k * DNZ;
					sum[field] += U[field][iii] * dV;
				}
			}
#pragma GCC ivdep
			for (integer k = HBW; k != HNX - HBW; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				const real rho = U[rho_i][iii];
				const real gx = G[gx_i][iii];
				const real gy = G[gy_i][iii];
				const real gz = G[gz_i][iii];
				const real x = X[XDIM][iii];
				const real y = X[YDIM][iii];
				const real z = X[ZDIM][iii];
				gx_sum += gx * rho * dV;
				gy_sum += gy * rho * dV;
				gz_sum += gz * rho * dV;
				glx_sum += (gy * z - gz * y) * rho * dV;
				gly_sum += (gx * z - gz * x) * rho * dV;
				glz_sum += (gx * y - gy * x) * rho * dV;
				lx_sum -= (z * U[sy_i][iii] - y * U[sz_i][iii]) * dV;
				ly_sum += (z * U[sx_i][iii] - x * U[sz_i][iii]) * dV;
				lz_sum -= (y * U[sx_i][iii] - x * U[sy_i][iii]) * dV;
				sx_sum += S[XDIM][iii] * dV;
				sy_sum += S[YDIM][iii] * dV;
				sz_sum += S[ZDIM][iii] * dV;
			}
		}
	}
	sum[egas_i] += HALF * sum[pot_i];

	FILE* fp = fopen("diag.dat", "at");
	fprintf(fp, "%18.10e ", double(t));
	for (integer field = 0; field != NF; ++field) {
		fprintf(fp, "%18.10e %18.10e ", double(sum[field]), double(U_out[field]));
	}
	fprintf(fp, "\n");
	fclose(fp);

	fp = fopen("grav.dat", "at");
	fprintf(fp, "%18.10e ", double(t));
	fprintf(fp, "%18.10e ", double(gx_sum));
	fprintf(fp, "%18.10e ", double(gy_sum));
	fprintf(fp, "%18.10e ", double(gz_sum));
	fprintf(fp, "%18.10e ", double(glx_sum));
	fprintf(fp, "%18.10e ", double(gly_sum));
	fprintf(fp, "%18.10e ", double(glz_sum));
	fprintf(fp, "\n");
	fclose(fp);

	fp = fopen("l.dat", "at");
	fprintf(fp, "%18.10e ", double(t));
	fprintf(fp, "%18.10e ", double(lx_sum));
	fprintf(fp, "%18.10e ", double(sx_sum));
	fprintf(fp, "%18.10e ", double(S_out[XDIM]));
	fprintf(fp, "%18.10e ", double(ly_sum));
	fprintf(fp, "%18.10e ", double(sy_sum));
	fprintf(fp, "%18.10e ", double(S_out[YDIM]));
	fprintf(fp, "%18.10e ", double(lz_sum));
	fprintf(fp, "%18.10e ", double(sz_sum));
	fprintf(fp, "%18.10e ", double(S_out[ZDIM]));
	fprintf(fp, "\n");
	fclose(fp);
}
