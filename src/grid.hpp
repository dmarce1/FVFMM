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
#include <list>
#include <set>


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <list>

enum gsolve_type {
	RHO, DRHODT
};

struct npair {
	integer lev;
	std::pair<integer, integer> loc;
};

typedef std::pair<integer, integer> dpair;

class grid {
public:
	typedef std::array<real, NDIM> xpoint;
	struct node_point {
		xpoint pt;
		integer index;
		bool operator==(const node_point& other) const;
		bool operator<(const node_point& other) const;
	};
private:
	static constexpr integer thread_cnt = 8;

	std::array<std::vector<real>, NF> U;
	std::vector<std::vector<multipole> > M;
	real dx, t;
	std::array<real, NDIM> xmin;
	integer step_num;
	integer nlevel;
	std::array<std::vector<real>, NF> U0;
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
	std::vector<std::vector<expansion> > L;
	std::vector<std::vector<expansion> > L_c;
	std::vector<npair> ilist_n;
	std::vector<dpair> ilist_d;
	std::array<std::vector<dpair>,NFACE> ilist_d_bnd;
	std::array<std::vector<npair>,NFACE> ilist_n_bnd;
	static bool float_eq(real a, real b);
	static bool xpoint_eq(const xpoint& a, const xpoint& b);
public:

	struct output_list_type {
		std::set<node_point> nodes;
		std::list<integer> zones;
		std::array<std::vector<real>,NF+NGF> data;
	};
	static void merge_output_lists(output_list_type& l1, const output_list_type& l2);

	real& hydro_value(integer, integer, integer, integer);
	real hydro_value(integer, integer, integer, integer) const;
	multipole& multipole_value(integer, integer, integer, integer);
	const multipole& multipole_value(integer, integer, integer, integer) const;
	bool refine_me(integer lev) const;
	integer level_count() const;
	void compute_ilist();
	void compute_dudt();
	void egas_to_etot();
	void etot_to_egas();
	void solve_gravity(gsolve_type = RHO);
	void compute_multipoles(gsolve_type);
	void compute_interactions(gsolve_type);
	void compute_boundary_interactions(gsolve_type, integer face);
	void compute_expansions(gsolve_type);
	real get_time() const;
	integer get_step() const;
	void save(const char* filename) const;
	void load(const char* filename);
	void diagnostics();
	grid(const std::function<std::vector<real>(real, real, real)>&, real dx = TWO / real(INX),
			std::array<real, NDIM> xmin = { -ONE, -ONE, -ONE });
	grid(real dx = TWO / real(INX), std::array<real, NDIM> xmin = { -ONE, -ONE, -ONE });
	void allocate();
	void reconstruct();
	void store();
	void restore();
	real compute_fluxes();
	void compute_sources();
	void boundaries();
	void set_physical_boundaries(integer);
	void next_u(integer rk, real dt);
	real step();
	static void output(const output_list_type&, const char*);
	output_list_type get_output_list() const;
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
};

#endif /* GRID_HPP_ */
