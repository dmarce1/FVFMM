/*
 * roe.cpp
 *
 *  Created on: May 28, 2015
 *      Author: dmarce1
 */

#include "roe.hpp"
#include <cmath>
#include <cassert>

const integer con_i = rho_i;
const integer acl_i = sx_i;
const integer acr_i = sy_i;
const integer sh1_i = sz_i;
const integer sh2_i = egas_i;

std::array<std::vector<real>, NF> roe_primitives(std::array<std::vector<real>, NF>& U) {
	const std::size_t sz = U[0].size();
	std::array<std::vector<real>, NF> V;
	for (integer f = 0; f != NF; ++f) {
		V[f].resize(sz);
	}
	for (std::size_t iii = 0; iii != sz; ++iii) {
		const real rho = U[rho_i][iii];
		const real vx = U[sx_i][iii] / rho;
		const real vy = U[sy_i][iii] / rho;
		const real vz = U[sz_i][iii] / rho;
		const real egas = U[egas_i][iii];
		const real tau = U[tau_i][iii];
		const real pot = U[pot_i][iii];
		const real ek = HALF * rho * (vx * vx + vy * vy + vz * vz);
		real ei = egas - ek;
		if (ei < de_switch2 * egas) {
			ei = std::pow(tau, fgamma);
		}
		assert( ei > real(0));
		assert( rho > real(0));
		const real p = (fgamma - ONE) * ei;
		const real h = (egas + p) / rho;
		V[rho_i][iii] = rho;
		V[vx_i][iii] = vx;
		V[vy_i][iii] = vy;
		V[vz_i][iii] = vz;
		V[h_i][iii] = h;
		V[tau_i][iii] = tau;
		V[pot_i][iii] = pot;
	}
	return V;
}

std::array<std::vector<real>, NF> roe_averages(const std::array<std::vector<real>, NF>& VL,
		const std::array<std::vector<real>, NF>& VR) {
	const std::size_t sz = VR[0].size();
	std::array<std::vector<real>, NF> V_;
	for (integer f = 0; f != NF; ++f) {
		V_[f].resize(sz);
	}
	for (std::size_t iii = 0; iii != sz; ++iii) {
		const real wr = std::sqrt(VR[rho_i][iii]);
		const real wl = std::sqrt(VL[rho_i][iii]);
		const real w0 = wr + wl;
		V_[rho_i][iii] = wr * wl;
		V_[vx_i][iii] = (wr * VR[vx_i][iii] + wl * VL[vx_i][iii]) / w0;
		V_[vy_i][iii] = (wr * VR[vy_i][iii] + wl * VL[vy_i][iii]) / w0;
		V_[vz_i][iii] = (wr * VR[vz_i][iii] + wl * VL[vz_i][iii]) / w0;
		V_[h_i][iii] = (wr * VR[h_i][iii] + wl * VL[h_i][iii]) / w0;
		V_[tau_i][iii] = (wr * VR[tau_i][iii] + wl * VL[tau_i][iii]) / w0;
		V_[pot_i][iii] = (wr * VR[pot_i][iii] + wl * VL[pot_i][iii]) / w0;
	}
	return V_;
}

real roe_fluxes(std::array<std::vector<real>, NF>& F, std::array<std::vector<real>, NF>& UL,
		std::array<std::vector<real>, NF>& UR, integer dimension) {
	const std::size_t sz = UL[0].size();

	/*auto phi0 = [](real lambda, real delta) {
	 if( std::abs(lambda) < delta) {
	 return (lambda*lambda + delta*delta)/(real(2.0)*delta);
	 } else {
	 return std::abs(lambda);
	 }
	 };*/

	const integer u_i = vx_i + dimension;
	const integer v_i = vx_i + (dimension == XDIM ? YDIM : XDIM);
	const integer w_i = vx_i + (dimension == ZDIM ? YDIM : ZDIM);
	/*
	 const auto VL = roe_primitives(UL);
	 const auto VR = roe_primitives(UR);
	 const auto V_ = roe_averages(VL, VR);
	 */
	real max_lambda = real(0);

	for (std::size_t iii = 0; iii != sz; ++iii) {

		/*	const real rho_r = VR[rho_i][iii];
		 const real u_r = VR[u_i][iii];
		 const real v_r = VR[v_i][iii];
		 const real w_r = VR[w_i][iii];
		 const real h_r = VR[h_i][iii];
		 const real tau_r = VR[tau_i][iii];
		 const real pot_r = VR[pot_i][iii];
		 const real p_r = std::max((fgamma - real(1)) / fgamma * rho_r * (h_r - HALF * (u_r * u_r + v_r * v_r + w_r * w_r)), 1.0e-5);
		 assert( rho_r > real(0));
		 assert( p_r > real(0));

		 const real rho_l = VL[rho_i][iii];
		 const real u_l = VL[u_i][iii];
		 const real v_l = VL[v_i][iii];
		 const real w_l = VL[w_i][iii];
		 const real h_l = VL[h_i][iii];
		 const real tau_l = VL[tau_i][iii];
		 const real pot_l = VL[pot_i][iii];
		 const real p_l = std::max((fgamma - real(1)) / fgamma * rho_l * (h_l - HALF * (u_l * u_l + v_l * v_l + w_l * w_l)), 1.0e-5);
		 assert( rho_l > real(0));
		 assert( p_l > real(0));

		 const real rho = V_[rho_i][iii];
		 const real u = V_[u_i][iii];
		 const real v = V_[v_i][iii];
		 const real w = V_[w_i][iii];
		 const real h = V_[h_i][iii];
		 const real p = std::max((fgamma - real(1)) / fgamma * rho * (h - HALF * (u * u + v * v + w * w)), 1.0e-5);
		 const real tau = V_[tau_i][iii];
		 const real pot = V_[pot_i][iii];
		 assert( rho > real(0));
		 const real c = std::sqrt(fgamma * p / rho);
		 assert( p > real(0));

		 const real drho = rho_r - rho_l;
		 const real du = u_r - u_l;
		 const real dv = v_r - v_l;
		 const real dw = w_r - w_l;
		 const real dp = p_r - p_l;
		 const real dtaus = tau_r / rho_r - tau_l / rho_l;
		 const real dpots = pot_r / rho_r - pot_l / rho_l;

		 const real lambda_con = phi0(u, du);
		 const real& lambda_sh1 = lambda_con;
		 const real& lambda_sh2 = lambda_con;
		 const real& lambda_tau = lambda_con;
		 const real& lambda_pot = lambda_con;
		 const real lambda_acl = phi0(u - c, du);
		 const real lambda_acr = phi0(u + c, du);

		 const real f_con = (drho - dp / (c * c)) * lambda_con;
		 const real f_acr = (du + dp / (rho * c)) * lambda_acr;
		 const real f_acl = (du - dp / (rho * c)) * lambda_acl;
		 const real f_sh1 = dv * lambda_sh1;
		 const real f_sh2 = dw * lambda_sh2;
		 const real f_taus = dtaus * lambda_tau;
		 const real f_pots = dpots * lambda_pot;

		 const real cinv = c;

		 const real f_roe_rho = f_con + HALF * rho * cinv * (f_acr - f_acl);
		 const real f_roe_u = u * f_con + HALF * rho * cinv * (f_acr * (u + c) - f_acl * (u - c));
		 const real f_roe_v = v * f_roe_rho + rho * f_sh1;
		 const real f_roe_w = w * f_roe_rho + rho * f_sh2;
		 const real f_roe_egas = f_con * HALF * (u * u + v * v + w * w) + rho * (v * f_sh1 + w * f_sh2)
		 + HALF * rho * cinv * (f_acr * (h + u * c) - f_acl * (h - u * c));
		 const real f_roe_tau = (tau / rho) * f_roe_rho + rho * f_taus;
		 const real f_roe_pot = (pot / rho) * f_roe_rho + rho * f_pots;

		 F[rho_i][iii] = HALF * (rho_r * u_r + rho_l * u_l - f_roe_rho);
		 F[u_i][iii] = HALF * (rho_r * u_r * u_r + p_r + rho_l * u_l * u_l + p_l - f_roe_u);
		 F[v_i][iii] = HALF * (rho_r * u_r * v_r + rho_l * u_l * v_l - f_roe_v);
		 F[w_i][iii] = HALF * (rho_r * u_r * w_r + rho_l * u_l * w_l - f_roe_w);
		 F[egas_i][iii] = HALF * (rho_r * u_r * h_r + rho_l * u_l * h_l - f_roe_egas);
		 F[tau_i][iii] = HALF * (tau_r * u_r + tau_l * u_l - f_roe_tau);
		 F[pot_i][iii] = HALF * (pot_r * u_r + pot_l * u_l - f_roe_pot);

		 const real this_max_lambda = std::max(std::max(lambda_acl, lambda_acr), lambda_con);*/

		const real v_r = UR[u_i][iii] / UR[rho_i][iii];
		real ei_r = UR[egas_i][iii]
				- HALF * (UR[u_i][iii] * UR[u_i][iii] + UR[v_i][iii] * UR[v_i][iii] + UR[w_i][iii] * UR[w_i][iii])
						/ UR[rho_i][iii];
		if (ei_r < de_switch2 * UR[egas_i][iii]) {
			ei_r = std::pow(UR[tau_i][iii], fgamma);
		}
		const real p_r = (fgamma - ONE) * ei_r;
		const real c_r = std::sqrt(fgamma * p_r / UR[rho_i][iii]);

		const real v_l = UL[u_i][iii] / UL[rho_i][iii];
		real ei_l = UL[egas_i][iii]
				- HALF * (UL[u_i][iii] * UL[u_i][iii] + UL[v_i][iii] * UL[v_i][iii] + UL[w_i][iii] * UL[w_i][iii])
						/ UL[rho_i][iii];
		if (ei_l < de_switch2 * UL[egas_i][iii]) {
			ei_l = std::pow(UL[tau_i][iii], fgamma);
		}
		const real p_l = (fgamma - ONE) * ei_l;
		const real c_l = std::sqrt(fgamma * p_l / UL[rho_i][iii]);

		const real a = std::max(std::abs(v_r) + c_r, std::abs(v_l) + c_l);

		for (integer field = 0; field != NF; ++field) {
			F[field][iii] = HALF
					* (UR[field][iii] * v_r + UL[field][iii] * v_l - a * (UR[field][iii] - UL[field][iii]));
		}
		F[u_i][iii] += HALF * (p_r + p_l);
		F[egas_i][iii] += HALF * (p_r * v_r + p_l * v_l);
		max_lambda = std::max(max_lambda, a);
	}

	return max_lambda;
}
