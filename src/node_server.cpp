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

integer node_server::local_node_count = 0;
hpx::lcos::local::spinlock node_server::timestep_lock;
integer node_server::timestep_node_count = 0;
real node_server::local_timestep = std::numeric_limits<real>::max();
std::shared_ptr<channel<real>> node_server::local_timestep_channel;
bool node_server::static_initialized(false);
std::atomic<integer> node_server::static_initializing(0);
std::list<const node_server*> node_server::local_node_list;
hpx::lcos::local::spinlock node_server::local_node_list_lock;

real node_server::step() {
	real dt = ZERO;

	std::vector<hpx::future<real>> child_futs(NCHILD);
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_futs[ci] = children[ci].step();
		}
	}
	compute_fmm(RHO);

	/*	real a;
	 const real dx = TWO / real(INX << my_location.level());
	 real cfl0 = cfl;
	 grid_ptr->store();

	 for (integer rk = 0; rk < NRK; ++rk) {
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
	 collect_hydro_boundaries(rk);
	 }

	 if (is_refined) {
	 hpx::wait_all(child_futs.begin(), child_futs.end());
	 }

	 ++step_num;*/
	return dt;
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
		(*i)->global_timestep_channel->set_value(t);
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
		for (integer rk = 0; rk != NRK; ++rk) {
			sibling_hydro_channels[rk][face] = std::make_shared<channel<std::vector<real>>>();
		}
		sibling_gravity_channels[face] = std::make_shared<channel<std::vector<real>>>();
	}
	for (integer ci = 0; ci != NCHILD; ++ci) {
		child_gravity_channels[ci] = std::make_shared<channel<multipole_pass_type>>();
	}
	parent_gravity_channel = std::make_shared<channel<expansion_pass_type>>();
	current_time = t;
	dx = TWO / real(INX << my_location.level());
	for (integer d = 0; d != NDIM; ++d) {
		xmin[d] = my_location.x_location(d);
	}
	const integer flags = ((my_location.level() == 0) ? GRID_IS_ROOT : 0) | GRID_IS_LEAF;
	if (current_time == ZERO) {
		grid_ptr = std::make_shared < grid > (problem, dx, xmin, flags);
	} else {
		grid_ptr = std::make_shared < grid > (dx, xmin, flags);
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
				children[ci] = hpx::invalid_id;
			}
			is_refined = false;
			const integer flags = ((my_location.level() == 0) ? GRID_IS_ROOT : 0) | GRID_IS_LEAF;
			grid_ptr = std::make_shared < grid > (dx, xmin, flags);
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
			const integer flags = (my_location.level() == 0) ? GRID_IS_ROOT : 0;
			grid_ptr = std::make_shared < grid > (dx, xmin, flags);

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
	printf("Setup done %llx %s\n", integer(my_location.unique_id()), my_location.to_str().c_str());
	for (integer si = 0; si != NFACE; ++si) {
		siblings[si] = sib_futs[si].get();
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

void node_server::compute_fmm(gsolve_type type) {
	multipole_pass_type m_in, m_out;
	expansion_pass_type l_out;
	m_out.first.resize(INX * INX * INX);
	m_out.second.resize(INX * INX * INX);
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			const integer x0 = ((ci >> 0) & 1) * INX / 2;
			const integer y0 = ((ci >> 0) & 1) * INX / 2;
			const integer z0 = ((ci >> 0) & 1) * INX / 2;
			m_in = child_gravity_channels[ci]->get();
			for (integer i = 0; i != INX / 2; ++i) {
				for (integer j = 0; j != INX / 2; ++j) {
					for (integer k = 0; k != INX / 2; ++k) {
						const integer ii = i * INX * INX / 4 + j * INX / 2 + k;
						const integer io = (i + x0) * INX * INX + (j + y0) * INX + k + z0;
						m_out.first[io] = m_in.first[ii];
						m_out.second[io] = m_in.second[ii];
					}
				}
			}
		}
		m_out = grid_ptr->compute_multipoles(type, &m_out);
	} else {
		m_out = grid_ptr->compute_multipoles(type);
	}
	hpx::future<void> parent_fut;
	std::vector<hpx::future<void>> child_futs(NCHILD);
	if (my_location.level() != 0) {
		parent.send_gravity_multipoles(std::move(m_out), my_location.get_child_index());
	} else {
		parent_fut = hpx::make_ready_future();
	}
	grid_ptr->compute_interactions(type);
	std::vector<hpx::future<void>> sib_futs(NFACE);
	for (integer si = 0; si != NFACE; ++si) {
		if (!my_location.is_physical_boundary(si)) {
			sib_futs[si] = siblings[si].send_gravity_boundary(get_gravity_boundary(si), si ^ 1);
		} else {
			sib_futs[si] = hpx::make_ready_future();
		}
	}
	for (integer si = 0; si != NFACE; ++si) {
		if (!my_location.is_physical_boundary(si)) {
			const std::vector<real> tmp = sibling_gravity_channels[si]->get();
			set_gravity_boundary(std::move(tmp), si);
		}
	}

	expansion_pass_type l_in;
	if (my_location.level() != 0) {
		l_in = parent_gravity_channel->get();
	}
	const expansion_pass_type ltmp = grid_ptr->compute_expansions(type, my_location.level() == 0 ? nullptr : &l_in);
	if (is_refined) {
		l_out.first.resize(INX * INX * INX / NCHILD);
		if (type == RHO) {
			l_out.second.resize(INX * INX * INX / NCHILD);
		}
		for (integer ci = 0; ci != NCHILD; ++ci) {
			const integer x0 = ((ci >> 0) & 1) * INX / 2;
			const integer y0 = ((ci >> 0) & 1) * INX / 2;
			const integer z0 = ((ci >> 0) & 1) * INX / 2;
			m_in = child_gravity_channels[ci]->get();
			for (integer i = 0; i != INX / 2; ++i) {
				for (integer j = 0; j != INX / 2; ++j) {
					for (integer k = 0; k != INX / 2; ++k) {
						const integer io = i * INX * INX / 4 + j * INX / 2 + k;
						const integer ii = (i + x0) * INX * INX + (j + y0) * INX + k + z0;
						l_out.first[io] = ltmp.first[ii];
						l_out.second[io] = ltmp.second[ii];
					}
				}
			}
			child_futs[ci] = children[ci].send_gravity_expansions(std::move(l_out));
		}
	} else {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_futs[ci] = hpx::make_ready_future();
		}
	}

	hpx::wait_all(sib_futs.begin(), sib_futs.end());
	hpx::wait_all(child_futs.begin(), child_futs.end());
	parent_fut.get();
}

