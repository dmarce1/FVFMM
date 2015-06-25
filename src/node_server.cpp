/*
 * node_server.cpp
 *
 *  Created on: Jun 11, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"
#include "problem.hpp"
#include <boost/thread/lock_guard.hpp>

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::simple_component<node_server>, node_server);
HPX_PLAIN_ACTION(node_server::get_local_timestep, get_local_timestep_action);
HPX_PLAIN_ACTION(node_server::set_global_timestep, set_global_timestep_action);
HPX_PLAIN_ACTION(node_server::output, output_action);

HPX_REGISTER_COMPONENT_MODULE();

integer node_server::local_node_count = 0;
hpx::lcos::local::spinlock node_server::timestep_lock;
integer node_server::timestep_node_count = 0;
real node_server::local_timestep = std::numeric_limits<real>::max();
std::shared_ptr<channel<real>> node_server::local_timestep_channel;
bool node_server::static_initialized(false);
std::atomic<integer> node_server::static_initializing(0);
std::list<const node_server*> node_server::local_node_list;
hpx::lcos::local::spinlock node_server::local_node_list_lock;

void node_server::collect_boundaries() {
	std::vector<hpx::future<void>> send_futs(NFACE);
	std::array<bool, NFACE> is_physical;
	for (integer face = 0; face != NFACE; ++face) {
		send_futs[face] = hpx::make_ready_future();
	}
	for (integer dim = 0; dim != NDIM; ++dim) {
		for (integer face = 2 * dim; face != 2 * dim + 2; ++face) {
			is_physical[face] = my_location.is_physical_boundary(face);
			if (!is_physical[face]) {
				const auto bdata = get_boundary(face);
				send_futs[face] = siblings[face].send_boundary(bdata, face ^ 1);
			}
		}
		for (integer face = 2 * dim; face != 2 * dim + 2; ++face) {
			if (!is_physical[face]) {
				const std::vector<real> bdata = sibling_channels[face]->get();
				set_boundary(bdata, face);
			} else {
				grid_ptr->set_physical_boundaries(face);
			}
		}
	}
	hpx::wait_all(send_futs.begin(), send_futs.end());
}

void node_server::step() {

	if (is_refined) {

		std::vector<hpx::future<void>> child_futs;
		child_futs.resize(NCHILD);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_futs[ci] = children[ci].step();
		}
		reduce_this_timestep(std::numeric_limits<real>::max());
		hpx::wait_all(child_futs.begin(), child_futs.end());

	} else {

		real dt, a;
		const real dx = TWO / real(INX << my_location.level());
		real cfl0 = cfl;
		grid_ptr->store();

		for (integer rk = 0; rk < 2; ++rk) {
			grid_ptr->reconstruct();
			a = grid_ptr->compute_fluxes();
			if (rk == 0) {
				dt = cfl0 * dx / a;
				reduce_this_timestep(dt);
			}
			grid_ptr->compute_sources();
			grid_ptr->compute_dudt();
			if (rk == 0) {
				dt = global_timestep_channel->get();
			}
			grid_ptr->next_u(0, dt);
			collect_boundaries();
		}
	}

	++step_num;

}

void node_server::reduce_this_timestep(double dt) {
	boost::lock_guard<hpx::lcos::local::spinlock> lock(timestep_lock);
	local_timestep = std::min(dt, local_timestep);
	++timestep_node_count;
	if (timestep_node_count == local_node_count) {
		local_timestep_channel->set_value(local_timestep);
		timestep_node_count = 0;
		if (hpx::get_locality_id() == 0) {
			hpx::thread([]() {
				const std::vector<hpx::id_type> ids = hpx::find_all_localities();
				const integer size = ids.size();
				std::vector<hpx::future<real>> futs(ids.size());
				for (integer i = 0; i != size; ++i) {
					futs[i] = hpx::async < get_local_timestep_action > (ids[i]);
				}
				real t = std::numeric_limits<real>::max();
				for (integer i = 0; i != size; ++i) {
					t = std::min(futs[i].get(), t);
				}
				for (integer i = 0; i != size; ++i) {
					hpx::apply<set_global_timestep_action>(ids[i],t);
				}
			}).detach();
		}
	}
}

real node_server::get_local_timestep() {
	return local_timestep_channel->get();
}

void node_server::set_global_timestep(real t) {
	boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
	for (auto i = local_node_list.begin(); i != local_node_list.end(); ++i) {
		if (!(*i)->is_refined) {
			(*i)->global_timestep_channel->set_value(t);
		}
	}
}

void node_server::static_initialize() {
	if (!static_initialized) {
		bool test = static_initializing++;
		if (!test) {
			local_timestep_channel = std::make_shared<channel<real>>();
			static_initialized = true;
		}
		while (!static_initialized) {
			hpx::this_thread::yield();
		}
	}
}

void node_server::initialize(real t) {
	step_num = 0;
	static_initialize();
	global_timestep_channel = std::make_shared<channel<real>>();
	{
		boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		++local_node_count;
		local_node_list.push_front(this);
		my_list_iterator = local_node_list.begin();
	}
	is_refined = false;
	siblings.resize(NFACE);
	for (integer face = 0; face != NFACE; ++face) {
		sibling_channels[face] = std::make_shared<channel<std::vector<real>>>();
	}
	current_time = t;
	const real dx = TWO / real(INX << my_location.level());
	std::array<real, NDIM> xmin;
	for (integer d = 0; d != NDIM; ++d) {
		xmin[d] = my_location.x_location(d);
	}
	if (current_time == ZERO) {
		grid_ptr = std::make_shared < grid > (problem, dx, xmin);
	} else {
		grid_ptr = std::make_shared < grid > (dx, xmin);
	}
}

node_server::~node_server() {
	{
		boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		--local_node_count;
		local_node_list.erase(my_list_iterator);
	}

}

node_server::node_server() {
	initialize(ZERO);
}

node_server::node_server(const node_location& loc, const node_client& parent_id, real t) :
		my_location(loc), parent(parent_id) {
	initialize(t);
}

integer node_server::regrid_gather() {
	integer count = integer(1);

	if (is_refined) {
		std::vector<hpx::future<integer>> futs(NCHILD);

		for (integer ci = 0; ci != NCHILD; ++ci) {
			futs[ci] = children[ci].regrid_gather();
		}

		for (integer ci = 0; ci != NCHILD; ++ci) {
			auto tmp = futs[ci].get();
			child_descendant_count[ci] = tmp;
			count += tmp;
		}

		if (count == 1) {
			std::vector<hpx::future<void>> futs(NCHILD);
			for (integer ci = 0; ci != NCHILD; ++ci) {
				futs[ci] = children[ci].unregister(my_location.get_child(ci));
			}
			is_refined = false;
			hpx::wait_all(futs.begin(), futs.end());
		}

	} else {
		if (grid_ptr->refine_me(my_location.level())) {
			count += NCHILD;

			/* Inefficient, only needs to be done once - rewrite*/
			me = my_location.get_client().get();

			std::vector < hpx::future < hpx::id_type >> idfuts(NCHILD);
			std::vector<hpx::future<void>> vfuts(NCHILD);
			children.resize(NCHILD);
			std::vector<node_location> clocs(NCHILD);
			for (integer ci = 0; ci != NCHILD; ++ci) {
				child_descendant_count[ci] = 1;
				clocs[ci] = my_location.get_child(ci);
				idfuts[ci] = hpx::new_ < node_server > (hpx::find_here(), clocs[ci], me, current_time);
			}
			for (integer ci = 0; ci != NCHILD; ++ci) {
				children[ci] = idfuts[ci].get();
				vfuts[ci] = children[ci].register_(clocs[ci]);
			}
			is_refined = true;
			hpx::wait_all(vfuts.begin(), vfuts.end());
		}
	}
	return count;
}

void node_server::regrid() {
	assert(grid_ptr!=nullptr);
	regrid_scatter(0, regrid_gather());
	assert(grid_ptr!=nullptr);
}

void node_server::regrid_scatter(integer a, integer total) {
	if (is_refined) {
		const auto localities = hpx::find_all_localities();
		++a;
		std::vector<hpx::future<void>> futs(NCHILD);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			const auto child_loc = localities[a * localities.size() / total];
			futs[ci] = children[ci].regrid_scatter(a, total);
			a += child_descendant_count[ci];
			//	hpx::components::migrate < node_server > (children[ci].get_gid(), child_loc);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
	std::vector<hpx::future<node_client>> sib_futs(NFACE);
	for (integer si = 0; si != NFACE; ++si) {
		if (!my_location.is_physical_boundary(si)) {
			sib_futs[si] = my_location.get_sibling(si).get_client();
		} else {
			sib_futs[si] = hpx::make_ready_future(node_client(hpx::id_type()));
		}
	}
	for (integer si = 0; si != NFACE; ++si) {
		siblings[si] = sib_futs[si].get();
	}
	printf("Done %llx\n", integer(my_location.unique_id()));
}

integer node_server::get_boundary_size(std::vector<std::array<integer, NDIM>>& lb,
		std::vector<std::array<integer, NDIM>>& ub, integer face, integer side) const {
	const integer lcnt = grid_ptr->level_count();
	integer hsize, msize, size, offset;
	lb.resize(lcnt);
	ub.resize(lcnt);
	size = 0;
	offset = (side == OUTER) ? HBW : 0;
	for (integer lev = 0; lev < lcnt; ++lev) {
		if (lev == 0) {
			msize = 1;
			hsize = NF;
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
			if (lev == 0) {
				hsize *= width;
			}
			msize *= width;
		}
		if (lev == 0) {
			size += hsize;
		}
		size += msize;
	}
	return size;
}

std::vector<real> node_server::get_boundary(integer face) {

	const integer lcnt = grid_ptr->level_count();
	std::vector<std::array<integer, NDIM>> lb, ub;
	std::vector<real> data;
	const integer size = get_boundary_size(lb, ub, face, INNER);
	data.resize(size);
	integer iter = 0;

	for (integer field = 0; field != NF; ++field) {
		for (integer i = lb[0][XDIM]; i < ub[0][XDIM]; ++i) {
			for (integer j = lb[0][YDIM]; j < ub[0][YDIM]; ++j) {
				for (integer k = lb[0][ZDIM]; k < ub[0][ZDIM]; ++k) {
					data[iter] = grid_ptr->hydro_value(field, i, j, k);
					++iter;
				}
			}
		}
	}

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

void node_server::set_boundary(const std::vector<real>& data, integer face) {
	const integer lcnt = grid_ptr->level_count();
	std::vector<std::array<integer, NDIM>> lb, ub;
	get_boundary_size(lb, ub, face, OUTER);
	integer iter = 0;

	for (integer field = 0; field != NF; ++field) {
		for (integer i = lb[0][XDIM]; i < ub[0][XDIM]; ++i) {
			for (integer j = lb[0][YDIM]; j < ub[0][YDIM]; ++j) {
				for (integer k = lb[0][ZDIM]; k < ub[0][ZDIM]; ++k) {
					grid_ptr->hydro_value(field, i, j, k) = data[iter];
					++iter;
				}
			}
		}
	}

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

void node_server::output(const std::string& filename) {
	const auto localities = hpx::find_all_localities();
	std::vector<hpx::future<void>> futs(localities.size());
	for (std::size_t i = 0; i != localities.size(); ++i) {
		if (localities[i] != hpx::find_here()) {
			futs[i] = hpx::async < output_action > (localities[i], filename);
		} else {
			futs[i] = hpx::make_ready_future();
		}
	}
	{
		grid::output_list_type olist;
		boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		for (auto i = local_node_list.begin(); i != local_node_list.end(); ++i) {
			if (!(*i)->is_refined) {
				grid::output_list_type this_list = (*i)->grid_ptr->get_output_list();
				grid::merge_output_lists(olist, std::move(this_list));
			}
		}
		grid::output(olist, filename.c_str());
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void node_server::recv_boundary(const std::vector<real>& bdata, integer face) {
	sibling_channels[face]->set_value(bdata);
}
