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

typedef node_server::regrid_gather_action regrid_gather_action_type;
typedef node_server::regrid_scatter_action regrid_scatter_action_type;
typedef node_server::send_hydro_boundary_action send_hydro_boundary_action_type;
typedef node_server::send_gravity_boundary_action send_gravity_boundary_action_type;
typedef node_server::send_gravity_multipoles_action send_gravity_multipoles_action_type;
typedef node_server::send_gravity_expansions_action send_gravity_expansions_action_type;
typedef node_server::step_action step_action_type;
typedef node_server::regrid_action regrid_action_type;
typedef node_server::solve_gravity_action solve_gravity_action_type;
typedef node_server::start_run_action start_run_action_type;
typedef node_server::copy_to_locality_action copy_to_locality_action_type;
typedef node_server::get_child_client_action get_child_client_action_type;
typedef node_server::form_tree_action form_tree_action_type;

HPX_REGISTER_ACTION (regrid_gather_action_type);
HPX_REGISTER_ACTION (regrid_scatter_action_type);
HPX_REGISTER_ACTION (send_hydro_boundary_action_type);
HPX_REGISTER_ACTION (send_gravity_boundary_action_type);
HPX_REGISTER_ACTION (send_gravity_multipoles_action_type);
HPX_REGISTER_ACTION (send_gravity_expansions_action_type);
HPX_REGISTER_ACTION (step_action_type);
HPX_REGISTER_ACTION (regrid_action_type);
HPX_REGISTER_ACTION (solve_gravity_action_type);
HPX_REGISTER_ACTION (start_run_action_type);
HPX_REGISTER_ACTION (copy_to_locality_action_type);
HPX_REGISTER_ACTION (get_child_client_action_type);
HPX_REGISTER_ACTION (form_tree_action_type);

integer node_server::local_node_count = 0;
hpx::lcos::local::spinlock node_server::timestep_lock;
integer node_server::timestep_node_count = 0;
real node_server::local_timestep = std::numeric_limits<real>::max();
std::shared_ptr<channel<real>> node_server::local_timestep_channel;
bool node_server::static_initialized(false);
std::atomic<integer> node_server::static_initializing(0);
std::list<const node_server*> node_server::local_node_list;
hpx::lcos::local::spinlock node_server::local_node_list_lock;

void node_server::form_tree(const hpx::id_type& self_gid, const hpx::id_type& parent_gid,
		const std::vector<hpx::shared_future<hpx::id_type>>& sib_gids) {
	for (integer si = 0; si != NFACE; ++si) {
		siblings[si] = sib_gids[si];
	}
	std::list<hpx::future<void>> cfuts;
	me = self_gid;
	parent = parent_gid;
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			std::vector < hpx::shared_future < hpx::id_type >> child_sibs(NFACE);
			for (integer d = 0; d != NDIM; ++d) {
				const integer flip = ci ^ (1 << d);
				const integer bit = (ci >> d) & 1;
				const integer other = 2 * d + bit;
				const integer thisf = 2 * d + (1 - bit);
				child_sibs[thisf] = hpx::make_ready_future(children[flip].get_gid());
				child_sibs[other] = siblings[other].get_child_client(flip);
			}
			cfuts.push_back(children[ci].form_tree(children[ci].get_gid(), me.get_gid(), child_sibs));
		}
		hpx::wait_all(cfuts.begin(), cfuts.end());
	}

}

hpx::id_type node_server::get_child_client(integer ci) {
	return children[ci].get_gid();
}

hpx::future<hpx::id_type> node_server::copy_to_locality(const hpx::id_type& id) {
	std::vector<hpx::id_type> cids;
	if (is_refined) {
		cids.resize(NCHILD);
		for (integer ci = 0; ci != NCHILD; ++ci) {
			cids[ci] = children[ci].get_gid();
		}
	}
	auto rc = hpx::new_ < node_server
			> (id, my_location, step_num, is_refined, current_time, child_descendant_count, *grid_ptr, cids);
	clear_family();
	return rc;
}

node_server::node_server(node_location&& _my_location, integer _step_num, bool _is_refined, real _current_time, std::array<integer,NCHILD>&& _child_d, grid&& _grid,const std::vector<hpx::id_type>& _c) {
	my_location = std::move(_my_location);
	initialize(_current_time);
	is_refined = _is_refined;
	step_num = _step_num;
	current_time = _current_time;
	grid_ptr = std::make_shared<grid>(std::move(_grid));
	if( is_refined ) {
		children.resize(NCHILD);
		for( integer ci = 0; ci != NCHILD; ++ci) {
			children[ci] = _c[ci];
		}
	}
	child_descendant_count = _child_d;
}

real node_server::step() {

	real dt = ZERO;

	std::list<hpx::future<real>> child_futs;
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_futs.push_back(children[ci].step());
		}
	}

	real a;
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
		compute_fmm(DRHODT, true, 2 * rk + 0);
		if (rk == 0) {
			dt = global_timestep_channel->get();
		}
		grid_ptr->next_u(rk, dt);
		compute_fmm(RHO, true, 2 * rk + 1);
		collect_hydro_boundaries(rk);
	}

	if (is_refined) {
		hpx::wait_all(child_futs.begin(), child_futs.end());
	}

	++step_num;
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
				std::list<hpx::future<real>> futs;
				for (integer i = 0; i != size; ++i) {
					futs.push_back(hpx::async < get_local_timestep_action > (ids[i]));
				}
				real t = std::numeric_limits<real>::max();
				for (auto i = futs.begin(); i != futs.end(); ++i) {
					t = std::min(i->get(), t);
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
			sibling_hydro_channels[rk][face] = std::make_shared<channel<std::vector<real>> >();
		}
	}
	for (integer c = 0; c != 4; ++c) {
		for (integer face = 0; face != NFACE; ++face) {
			sibling_gravity_channels[c][face] = std::make_shared<channel<std::vector<real>> >();
		}
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_gravity_channels[c][ci] = std::make_shared<channel<multipole_pass_type>>();
		}
		parent_gravity_channel[c] = std::make_shared<channel<expansion_pass_type>>();
	}
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
	printf("Creating grid at %i: %i %i %i w on locality %i\n", int(my_location.level()), int(my_location[0]),
			int(my_location[1]), int(my_location[2]), int(hpx::get_locality_id()));
}

node_server::~node_server() {

	boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
	--local_node_count;
	local_node_list.erase(my_list_iterator);
	printf("Destroying grid at %i: %i %i %i w on locality %i\n", int(my_location.level()), int(my_location[0]),
			int(my_location[1]), int(my_location[2]), int(hpx::get_locality_id()));
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
		std::vector<hpx::shared_future<integer>> futs(NCHILD);

		for (integer ci = 0; ci != NCHILD; ++ci) {
			futs[ci] = children[ci].regrid_gather();
		}

		for (integer ci = 0; ci != NCHILD; ++ci) {
			auto tmp = futs[ci].get();
			child_descendant_count[ci] = tmp;
			count += tmp;
		}

		if (count == 1) {
			for (integer ci = 0; ci != NCHILD; ++ci) {
				children[ci] = node_client();
			}
			is_refined = false;
			const integer flags = ((my_location.level() == 0) ? GRID_IS_ROOT : 0) | GRID_IS_LEAF;
			grid_ptr = std::make_shared < grid > (dx, xmin, flags);
		}

	} else {
		if (grid_ptr->refine_me(my_location.level())) {
			count += NCHILD;

			/* Inefficient, only needs to be done once - rewrite*/
			me = my_location.get_client().get();

			children.resize(NCHILD);
			std::vector<node_location> clocs(NCHILD);
			for (integer ci = 0; ci != NCHILD; ++ci) {
				child_descendant_count[ci] = 1;
				clocs[ci] = my_location.get_child(ci);
				children[ci] = node_client::create(hpx::find_here(), clocs[ci], me, current_time);
			}
			is_refined = true;
			const integer flags = (my_location.level() == 0) ? GRID_IS_ROOT : 0;
			grid_ptr = std::make_shared < grid > (dx, xmin, flags);
		}
	}
	return count;
}

void node_server::regrid() {
	assert(grid_ptr!=nullptr);
	regrid_scatter(0, regrid_gather());
	assert(grid_ptr!=nullptr);
}

void node_server::regrid_scatter(integer a_, integer total) {
	std::list<hpx::future<void>> futs;
//	std::vector<hpx::future<void>> futs2(NCHILD);
	if (is_refined) {
		integer a = a_;
		const auto localities = hpx::find_all_localities();
		++a;
		for (integer ci = 0; ci != NCHILD; ++ci) {
			const integer loc_index = a * localities.size() / total;
			const auto child_loc = localities[loc_index];
			a += child_descendant_count[ci];
			children[ci] = children[ci].copy_to_locality(child_loc);
		}
		a = a_ + 1;
		for (integer ci = 0; ci != NCHILD; ++ci) {
			futs.push_back(children[ci].regrid_scatter(a, total));
			a += child_descendant_count[ci];
		}
	}
	clear_family();
	if (is_refined) {
		hpx::wait_all(futs.begin(), futs.end());
	}
}

void node_server::clear_family() {
	for (integer si = 0; si != NFACE; ++si) {
		siblings[si] = hpx::invalid_id;
	}
	parent = hpx::invalid_id;
	me = hpx::invalid_id;
}

grid::output_list_type node_server::output(const std::string& filename) {
	const auto localities = hpx::find_all_localities();
	std::list < hpx::future < grid::output_list_type >> futs;
	if (hpx::get_locality_id() == 0) {
		for (std::size_t i = 0; i != localities.size(); ++i) {
			if (localities[i] != hpx::find_here()) {
				futs.push_back(hpx::async < output_action > (localities[i], filename));
			}
		}
	}
	grid::output_list_type olist;
	{
		boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		for (auto i = local_node_list.begin(); i != local_node_list.end(); ++i) {
			if (!(*i)->is_refined) {
				grid::output_list_type this_list = (*i)->grid_ptr->get_output_list();
				grid::merge_output_lists(olist, std::move(this_list));
			}
		}
	}
	if (hpx::get_locality_id() == 0) {
		for (auto i = futs.begin(); i != futs.end(); ++i) {
			grid::merge_output_lists(olist, i->get());
		}
		grid::output(olist, filename.c_str());
	}
	return olist;
}

void node_server::solve_gravity(bool ene, integer c) {
//	printf("%s\n", my_location.to_str().c_str());

	std::list<hpx::future<void>> child_futs;
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_futs.push_back(children[ci].solve_gravity(ene, c));
		}
	}
	compute_fmm(RHO, ene, c);
	if (is_refined) {
		hpx::wait_all(child_futs.begin(), child_futs.end());
	}
}

void node_server::compute_fmm(gsolve_type type, bool energy_account, integer gchannel) {
	if (energy_account) {
		grid_ptr->egas_to_etot();
	}
	multipole_pass_type m_in, m_out;
	expansion_pass_type l_out;
	m_out.first.resize(INX * INX * INX);
	m_out.second.resize(INX * INX * INX);
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			const integer x0 = ((ci >> 0) & 1) * INX / 2;
			const integer y0 = ((ci >> 1) & 1) * INX / 2;
			const integer z0 = ((ci >> 2) & 1) * INX / 2;
			m_in = child_gravity_channels[gchannel][ci]->get();
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
	std::list<hpx::future<void>> child_futs;
	if (my_location.level() != 0) {
		parent_fut = parent.send_gravity_multipoles(std::move(m_out), my_location.get_child_index(), gchannel);
	} else {
		parent_fut = hpx::make_ready_future();
	}
	grid_ptr->compute_interactions(type);

	for (integer dim = 0; dim != NDIM; ++dim) {
		std::list<hpx::future<void>> sib_futs;
		for (integer si = 2 * dim; si != 2 * (dim + 1); ++si) {
			if (!my_location.is_physical_boundary(si)) {
				sib_futs.push_back(siblings[si].send_gravity_boundary(get_gravity_boundary(si), si ^ 1, gchannel));
			} else {
				sib_futs.push_back(hpx::make_ready_future());
			}
		}
		for (integer si = 2 * dim; si != 2 * (dim + 1); ++si) {
			if (!my_location.is_physical_boundary(si)) {
				const std::vector<real> tmp = sibling_gravity_channels[gchannel][si]->get();
				set_gravity_boundary(std::move(tmp), si);
				grid_ptr->compute_boundary_interactions(type, si);
			}
		}
		for (auto si = sib_futs.begin(); si != sib_futs.end(); ++si) {
			si->get();
		}
	}

	expansion_pass_type l_in;
	if (my_location.level() != 0) {
		l_in = parent_gravity_channel[gchannel]->get();
	}
	const expansion_pass_type ltmp = grid_ptr->compute_expansions(type, my_location.level() == 0 ? nullptr : &l_in);
	if (is_refined) {
		l_out.first.resize(INX * INX * INX / NCHILD);
		if (type == RHO) {
			l_out.second.resize(INX * INX * INX / NCHILD);
		}
		for (integer ci = 0; ci != NCHILD; ++ci) {
			const integer x0 = ((ci >> 0) & 1) * INX / 2;
			const integer y0 = ((ci >> 1) & 1) * INX / 2;
			const integer z0 = ((ci >> 2) & 1) * INX / 2;
			for (integer i = 0; i != INX / 2; ++i) {
				for (integer j = 0; j != INX / 2; ++j) {
					for (integer k = 0; k != INX / 2; ++k) {
						const integer io = i * INX * INX / 4 + j * INX / 2 + k;
						const integer ii = (i + x0) * INX * INX + (j + y0) * INX + k + z0;
						l_out.first[io] = ltmp.first[ii];
						if (type == RHO) {
							l_out.second[io] = ltmp.second[ii];
						}
					}
				}
			}
			child_futs.push_back(children[ci].send_gravity_expansions(std::move(l_out), gchannel));
		}
	}
	parent_fut.get();
	hpx::wait_all(child_futs.begin(), child_futs.end());

	if (energy_account) {
		grid_ptr->etot_to_egas();
	}
}

