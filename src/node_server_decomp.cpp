/*
 * node_server_decomp.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"

void node_server::collect_hydro_boundaries(integer rk) {
	std::vector<hpx::future<void>> send_futs(NFACE);
	std::array<bool, NFACE> is_physical;
	for (integer face = 0; face != NFACE; ++face) {
		send_futs[face] = hpx::make_ready_future();
	}
	for (integer dim = 0; dim != NDIM; ++dim) {
		for (integer face = 2 * dim; face != 2 * dim + 2; ++face) {
			is_physical[face] = my_location.is_physical_boundary(face);
			if (!is_physical[face]) {
				const auto bdata = get_hydro_boundary(face);
				send_futs[face] = siblings[face].send_hydro_boundary(bdata, rk, face ^ 1);
			}
		}
		for (integer face = 2 * dim; face != 2 * dim + 2; ++face) {
			if (!is_physical[face]) {
				const std::vector<real> bdata = sibling_hydro_channels[rk][face]->get();
				set_hydro_boundary(bdata, face);
			} else {
				grid_ptr->set_physical_boundaries(face);
			}
		}
	}
	hpx::wait_all(send_futs.begin(), send_futs.end());
}

integer node_server::get_hydro_boundary_size(std::array<integer, NDIM>& lb, std::array<integer, NDIM>& ub, integer face,
		integer side) const {
	integer hsize, size, offset;
	size = 0;
	offset = (side == OUTER) ? HBW : 0;
	hsize = NF;
	for (integer d = 0; d != NDIM; ++d) {
		const integer nx = INX + 2 * HBW;
		if (d < face / 2) {
			lb[d] = 0;
			ub[d] = nx;
		} else if (d > face / 2) {
			lb[d] = HBW;
			ub[d] = nx - HBW;
		} else if (face % 2 == 0) {
			lb[d] = HBW - offset;
			ub[d] = 2 * HBW - offset;
		} else {
			lb[d] = nx - 2 * HBW + offset;
			ub[d] = nx - HBW + offset;
		}
		const integer width = ub[d] - lb[d];
		hsize *= width;
	}
	size += hsize;
	return size;
}

integer node_server::get_gravity_boundary_size(std::vector<std::array<integer, NDIM>>& lb,
		std::vector<std::array<integer, NDIM>>& ub, integer face, integer side) const {
	const integer lcnt = grid_ptr->level_count();
	integer msize, size, offset;
	lb.resize(lcnt);
	ub.resize(lcnt);
	size = 0;
	offset = (side == OUTER) ? HBW : 0;
	for (integer lev = 0; lev < lcnt; ++lev) {
		if (lev == 0) {
			msize = 1;
		} else {
			msize = 20;
		}
		for (integer d = 0; d != NDIM; ++d) {
			const integer nx = (INX >> lev) + 2 * HBW;
			if (d < face / 2) {
				lb[lev][d] = 0;
				ub[lev][d] = nx;
			} else if (d > face / 2) {
				lb[lev][d] = HBW;
				ub[lev][d] = nx - HBW;
			} else if (face % 2 == 0) {
				lb[lev][d] = HBW - offset;
				ub[lev][d] = 2 * HBW - offset;
			} else {
				lb[lev][d] = nx - 2 * HBW + offset;
				ub[lev][d] = nx - HBW + offset;
			}
			const integer width = ub[lev][d] - lb[lev][d];
			msize *= width;
		}
		size += msize;
	}
	return size;
}

std::vector<real> node_server::get_hydro_boundary(integer face) {

	std::array<integer, NDIM> lb, ub;
	std::vector<real> data;
	const integer size = get_hydro_boundary_size(lb, ub, face, INNER);
	data.resize(size);
	integer iter = 0;

	for (integer field = 0; field != NF; ++field) {
		for (integer i = lb[XDIM]; i < ub[XDIM]; ++i) {
			for (integer j = lb[YDIM]; j < ub[YDIM]; ++j) {
				for (integer k = lb[ZDIM]; k < ub[ZDIM]; ++k) {
					data[iter] = grid_ptr->hydro_value(field, i, j, k);
					++iter;
				}
			}
		}
	}

	return data;
}

std::vector<real> node_server::get_gravity_boundary(integer face) {

	const integer lcnt = grid_ptr->level_count();
	std::vector<std::array<integer, NDIM>> lb, ub;
	std::vector<real> data;
	const integer size = get_gravity_boundary_size(lb, ub, face, INNER);
	data.resize(size);
	integer iter = 0;

	for (integer lev = 0; lev < lcnt; ++lev) {
		for (integer i = lb[lev][XDIM]; i < ub[lev][XDIM]; ++i) {
			for (integer j = lb[lev][YDIM]; j < ub[lev][YDIM]; ++j) {
				for (integer k = lb[lev][ZDIM]; k < ub[lev][ZDIM]; ++k) {
					const auto& m = grid_ptr->multipole_value(lev, i, j, k);
					const integer top = lev == 0 ? 1 : 20;
					for (integer l = 0; l < top; ++l) {
						data[iter] = m.ptr()[l];
						++iter;
					}
				}
			}
		}
	}

	return data;
}

void node_server::set_hydro_boundary(const std::vector<real>& data, integer face) {
	std::array<integer, NDIM> lb, ub;
	get_hydro_boundary_size(lb, ub, face, OUTER);
	integer iter = 0;

	for (integer field = 0; field != NF; ++field) {
		for (integer i = lb[XDIM]; i < ub[XDIM]; ++i) {
			for (integer j = lb[YDIM]; j < ub[YDIM]; ++j) {
				for (integer k = lb[ZDIM]; k < ub[ZDIM]; ++k) {
					grid_ptr->hydro_value(field, i, j, k) = data[iter];
					++iter;
				}
			}
		}
	}
}

void node_server::set_gravity_boundary(const std::vector<real>& data, integer face) {
	const integer lcnt = grid_ptr->level_count();
	std::vector<std::array<integer, NDIM>> lb, ub;
	get_gravity_boundary_size(lb, ub, face, OUTER);
	integer iter = 0;

	for (integer lev = 0; lev < lcnt; ++lev) {
		for (integer i = lb[lev][XDIM]; i < ub[lev][XDIM]; ++i) {
			for (integer j = lb[lev][YDIM]; j < ub[lev][YDIM]; ++j) {
				for (integer k = lb[lev][ZDIM]; k < ub[lev][ZDIM]; ++k) {
					auto& m = grid_ptr->multipole_value(lev, i, j, k);
					const integer top = lev == 0 ? 1 : 20;
					for (integer l = 0; l < top; ++l) {
						m.ptr()[l] = data[iter];
						++iter;
					}
				}
			}
		}
	}
}

void node_server::recv_gravity_multipoles(multipole_pass_type&& v, integer ci) {
	child_gravity_channels[ci]->set_value(std::move(v));
}

void node_server::recv_gravity_expansions(expansion_pass_type&& v) {
	parent_gravity_channel->set_value(std::move(v));
}

void node_server::recv_hydro_boundary(std::vector<real>&& bdata, integer rk, integer face) {
	sibling_hydro_channels[rk][face]->set_value(std::move(bdata));
}

void node_server::recv_gravity_boundary(std::vector<real>&& bdata, integer face) {
	sibling_gravity_channels[face]->set_value(std::move(bdata));
}
