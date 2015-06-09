/*
 * grid.hpp
 *
 *  Created on: May 26, 2015
 *      Author: dmarce1
 */

#ifndef GRID_HPP_
#define GRID_HPP_

#include "defs.hpp"
#include "roe.hpp"
#include "taylor.hpp"
#include "space_vector.hpp"
#include <functional>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <list>

enum gsolve_type {RHO, DRHODT};

struct npair {
	integer lev;
	std::pair<integer,integer> loc;
};

typedef std::pair<integer,integer> dpair;

class grid {
private:
	real dx, t;
	integer step_num;
	integer nlevel;
	std::array<std::vector<real>, NF> U0;
	std::array<std::vector<real>, NF> U;
	std::array<std::vector<real>, NDIM> S0;
	std::array<std::vector<real>, NDIM> S;
	std::array<std::vector<real>, NF> src;
	std::vector<real> U_out;
	std::vector<real> U_out0;
	std::vector<real> S_out;
	std::vector<real> S_out0;
	std::array<std::vector<real>, NF> dUdt;
	std::array<std::array<std::vector<real>, NF>, NFACE> Uf;
	std::array<std::array<std::vector<real>, NF>, NDIM> F;
	std::array<std::vector<real>, NDIM> X;
	std::array<std::vector<real>, NGF> G;
	std::array<std::vector<real>, NGF> G0;
	std::vector<real> dphi_dt;
	std::vector<std::vector<space_vector> > com;
	std::vector<std::vector<multipole> > M;
	std::array<std::vector<std::vector<expansion> >, NGF> L;
	std::list<npair> ilist_n;
	std::list<dpair> ilist_d;
public:
	void compute_ilist();
	template<class Archive>
	void serialize(Archive& arc, const unsigned) {
		arc & U;
		arc & S;
		arc & G;
		arc & U_out;
		arc & S_out;
		arc & dx;
		arc & t;
		arc & step_num;
	}
	void compute_dudt();
	void egas_to_etot();
	void etot_to_egas();
	void solve_gravity(gsolve_type=RHO);
	void compute_multipoles(gsolve_type);
	void compute_interactions(gsolve_type);
	void compute_expansions(gsolve_type);
	real get_time() const;
	integer get_step() const;
	void save(const char* filename) const;
	void load(const char* filename);
	void diagnostics();
	grid(const std::function<std::vector<real>(real, real, real)>&);
	void reconstruct();
	void store();
	void restore();
	real compute_fluxes();
	void compute_sources();
	void boundaries();
	void next_u(integer rk, real dt);
	real step();
	void output(const char*);
};

#endif /* GRID_HPP_ */
