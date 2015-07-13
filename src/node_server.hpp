/*
 * node_server.hpp
 *
 *  Created on: Jun 11, 2015
 *      Author: dmarce1
 */

#ifndef NODE_SERVER_HPP_
#define NODE_SERVER_HPP_

#include "defs.hpp"
#include "node_location.hpp"
#include "node_client.hpp"
#include "grid.hpp"
#include "channel.hpp"

const integer INNER = 0;
const integer OUTER = 1;

class node_server: public hpx::components::simple_component_base<node_server> {
private:
	node_location my_location;
	integer step_num;
	real current_time;
	std::shared_ptr<grid> grid_ptr;
	bool is_refined;
	std::array<integer, NVERTEX> child_descendant_count;
	std::array<real, NDIM> xmin;
	real dx;
	node_client me;
	node_client parent;
	std::vector<node_client> siblings;
	std::vector<node_client> children;
public:
	template<class Archive>
	void serialize(Archive& arc, unsigned) {
		arc & my_location;
		arc & step_num;
		arc & is_refined;
		arc & current_time;
		arc & xmin;
		arc & dx;
		arc & *grid_ptr;
	}

	node_server(node_location&&, integer, bool, real, std::array<integer,NCHILD>&&, grid&&, const std::vector<hpx::id_type>&);
private:
	std::array<std::array<std::shared_ptr<channel<std::vector<real>>> ,NFACE>,NRK> sibling_hydro_channels;
	std::shared_ptr<channel<expansion_pass_type>> parent_gravity_channel;
	std::array<std::shared_ptr<channel<std::vector<real>>> ,NFACE> sibling_gravity_channels;
	std::array<std::shared_ptr<channel<multipole_pass_type>>, NCHILD> child_gravity_channels;
	std::shared_ptr<channel<real>> global_timestep_channel;

	std::list<const node_server*>::iterator my_list_iterator;

	static std::list<const node_server*> local_node_list;
	static hpx::lcos::local::spinlock local_node_list_lock;

	static integer local_node_count;
	static hpx::lcos::local::spinlock timestep_lock;
	static integer timestep_node_count;
	static real local_timestep;
	static bool static_initialized;
	static std::atomic<integer> static_initializing;
	static std::shared_ptr<channel<real>> local_timestep_channel;

	static void reduce_this_timestep(double dt);

	void initialize(real);
	void collect_hydro_boundaries(integer rk);
	static void static_initialize();
	void clear_family();

public:

	static grid::output_list_type output(const std::string&);
	static real get_local_timestep();
	static void set_global_timestep(real);

	node_server();
	~node_server();
	node_server( const node_server& other);
	node_server( node_server&& other);
	node_server(const node_location&, const node_client& parent_id, real);
	integer get_boundary_size(std::array<integer, NDIM>&, std::array<integer, NDIM>&, integer,
			integer) const;
	void set_hydro_boundary(const std::vector<real>&, integer face);
	std::vector<real> get_hydro_boundary(integer face);

	void set_gravity_boundary(const std::vector<real>&, integer face);
	std::vector<real> get_gravity_boundary(integer face);

	integer regrid_gather();
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid_gather, regrid_gather_action);

	void regrid_scatter(integer, integer);
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid_scatter, regrid_scatter_action);

	void recv_hydro_boundary( std::vector<real>&&, integer rk, integer face);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_hydro_boundary, send_hydro_boundary_action);

	void recv_gravity_boundary( std::vector<real>&&, integer face);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_gravity_boundary, send_gravity_boundary_action);

	void recv_gravity_multipoles( multipole_pass_type&&, integer ci);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_gravity_multipoles, send_gravity_multipoles_action);

	void recv_gravity_expansions(expansion_pass_type&&);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_gravity_expansions, send_gravity_expansions_action);

	real step();
	HPX_DEFINE_COMPONENT_ACTION(node_server, step, step_action);

	void regrid();
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid, regrid_action);

	void compute_fmm(gsolve_type gs, bool energy_account = true);

	void solve_gravity(bool ene=true);
	HPX_DEFINE_COMPONENT_ACTION(node_server, solve_gravity, solve_gravity_action);

	void start_run();
	HPX_DEFINE_COMPONENT_ACTION(node_server, start_run, start_run_action);

	hpx::future<hpx::id_type> copy_to_locality(const hpx::id_type& );
	HPX_DEFINE_COMPONENT_ACTION(node_server, copy_to_locality, copy_to_locality_action);

//	void find_family();
//	HPX_DEFINE_COMPONENT_ACTION(node_server, find_family, find_family_action);

	hpx::id_type get_child_client(integer ci);
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_child_client, get_child_client_action);

	std::vector<hpx::shared_future<hpx::id_type>> get_sibling_clients( integer ci);
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_sibling_clients, get_sibling_clients_action);

	void form_tree(const hpx::id_type&, const hpx::id_type&, const std::vector<hpx::shared_future<hpx::id_type>>& );
	HPX_DEFINE_COMPONENT_ACTION(node_server, form_tree, form_tree_action);

};

//HPX_REGISTER_ACTION_DECLARATION (node_server::start_run_action);
//HPX_REGISTER_ACTION_DECLARATION (node_server::step_action);
//HPX_REGISTER_ACTION_DECLARATION (node_server::solve_gravity_action);

#endif /* NODE_SERVER_HPP_ */
