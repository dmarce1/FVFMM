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
#include <atomic>

const integer INNER = 0;
const integer OUTER = 1;

struct diagnostics_t {
	std::vector<real> grid_sum;
	std::vector<real> outflow_sum;
	std::vector<real> l_sum;
	std::vector<real> s_sum;
	diagnostics_t() : grid_sum(NF,ZERO), outflow_sum(NF,ZERO), l_sum(NDIM), s_sum(NDIM){}
	diagnostics_t& operator+=(const diagnostics_t& other) {
		for( integer f = 0; f != NF; ++f) {
			grid_sum[f] += other.grid_sum[f];
			outflow_sum[f] += other.outflow_sum[f];
		}
		for( integer d = 0; d != NDIM; ++d) {
			l_sum[d] += other.l_sum[d];
			s_sum[d] += other.s_sum[d];
		}
		return *this;
	}
	diagnostics_t& operator=(const diagnostics_t& other) {
			grid_sum = other.grid_sum;
			outflow_sum = other.outflow_sum;
			s_sum = other.s_sum;
			l_sum = other.l_sum;
			return *this;
	};
	template<class Arc>
	void serialize( Arc& arc, const unsigned ) {
		arc & grid_sum;
		arc & outflow_sum;
		arc & l_sum;
		arc & s_sum;
	}
};

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
	real get_time() const {
		return current_time;
	}
	node_server& operator=( node_server&& );

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



	static std::pair<integer,std::size_t> save(const std::string& filename);

	static void load(const std::string& filename, node_client root);

	node_server(node_location&&, integer, bool, real, std::array<integer,NCHILD>&&, grid&&, const std::vector<hpx::id_type>&);
	node_server(node_server&& other);

	std::size_t load_me( FILE* fp, integer&);
	std::size_t save_me( FILE* fp ) const;
private:
	std::array<std::array<std::shared_ptr<channel<std::vector<real>>> ,NFACE>,NRK> sibling_hydro_channels;
	std::array<std::shared_ptr<channel<expansion_pass_type>>,4> parent_gravity_channel;
	std::array<std::array<std::shared_ptr<channel<std::vector<real>>> ,NFACE>,4> sibling_gravity_channels;
	std::array<std::array<std::shared_ptr<channel<std::vector<real> >>, NCHILD>,NRK> child_hydro_channels;
	std::array<std::array<std::shared_ptr<channel<multipole_pass_type>>, NCHILD>,4> child_gravity_channels;
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
	hpx::future<void> exchange_interlevel_hydro_data(integer rk);
	static void static_initialize();
	void clear_family();

public:

	static void output_form();

	static grid::output_list_type output_collect(const std::string&);

	static 	std::list<grid::output_list_type> olists;


	static real get_local_timestep();
	static void set_global_timestep(real);

	node_server();
	~node_server();
	node_server( const node_server& other);
	node_server(const node_location&, const node_client& parent_id, real);
	std::vector<real> restricted_grid() const;
	void load_from_restricted_child(const std::vector<real>&, integer);
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

	void recv_hydro_children( std::vector<real>&&, integer rk, integer ci);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_hydro_children, send_hydro_children_action);

	void recv_gravity_boundary( std::vector<real>&&, integer face, integer);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_gravity_boundary, send_gravity_boundary_action);

	void recv_gravity_multipoles( multipole_pass_type&&, integer ci, integer);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_gravity_multipoles, send_gravity_multipoles_action);

	void recv_gravity_expansions(expansion_pass_type&&, integer);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_gravity_expansions, send_gravity_expansions_action);

	real step();
	HPX_DEFINE_COMPONENT_ACTION(node_server, step, step_action);

	void regrid();
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid, regrid_action);

	void compute_fmm(gsolve_type gs, bool energy_account, integer c);

	void solve_gravity(bool ene, integer c);
	HPX_DEFINE_COMPONENT_ACTION(node_server, solve_gravity, solve_gravity_action);

	void start_run();
	HPX_DEFINE_COMPONENT_ACTION(node_server, start_run, start_run_action);

	hpx::future<hpx::id_type> copy_to_locality(const hpx::id_type& );
	HPX_DEFINE_COMPONENT_ACTION(node_server, copy_to_locality, copy_to_locality_action);

	hpx::id_type get_child_client(integer ci);
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_child_client, get_child_client_action);

	std::vector<hpx::shared_future<hpx::id_type>> get_sibling_clients( integer ci);
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_sibling_clients, get_sibling_clients_action);


	void form_tree(const hpx::id_type&, const hpx::id_type&, const std::vector<hpx::shared_future<hpx::id_type>>& );
	HPX_DEFINE_COMPONENT_ACTION(node_server, form_tree, form_tree_action);

	hpx::id_type load_node( std::size_t, const std::string&, const node_location&, const hpx::id_type& );
	HPX_DEFINE_COMPONENT_ACTION(node_server, load_node, load_node_action);

	node_server* get_ptr();
	HPX_DEFINE_COMPONENT_ACTION(node_server, get_ptr, get_ptr_action);

	diagnostics_t diagnostics() const;
	HPX_DEFINE_COMPONENT_ACTION(node_server, diagnostics, diagnostics_action);

};


HPX_REGISTER_ACTION_DECLARATION( node_server::regrid_gather_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::regrid_scatter_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::send_hydro_boundary_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::send_gravity_boundary_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::send_gravity_multipoles_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::send_gravity_expansions_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::send_hydro_children_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::step_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::regrid_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::solve_gravity_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::start_run_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::copy_to_locality_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::get_child_client_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::form_tree_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::get_ptr_action);
HPX_REGISTER_ACTION_DECLARATION( node_server::diagnostics_action);

HPX_DEFINE_PLAIN_ACTION(node_server::save, save_action);

#endif /* NODE_SERVER_HPP_ */
