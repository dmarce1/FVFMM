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
	integer step_num;
	real current_time;
	std::shared_ptr<grid> grid_ptr;
	node_location my_location;
	node_client me;
	node_client parent;
	std::vector<node_client> siblings;
	std::vector<node_client> children;
	bool is_refined;
	std::array<integer, NVERTEX> child_descendant_count;
	integer locality_id;
	std::array<std::array<std::shared_ptr<channel<std::vector<real>>> ,NFACE>,NRK> sibling_channels;
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
	std::shared_ptr<channel<real>> global_timestep_channel;

	static void reduce_this_timestep(double dt);

	void initialize(real);
	void collect_boundaries(integer rk);
	static void static_initialize();
public:

	static void output(const std::string&);
	static real get_local_timestep();
	static void set_global_timestep(real);

	template<class Archive>
	void serialize(Archive& arc, unsigned) {
		arc & my_location;
		arc & step_num;
		arc & me;
		arc & parent;
		arc & siblings;
		arc & children;
		arc & is_refined;
		arc & child_descendant_count;
		arc & locality_id;
		arc & *grid_ptr;
	}

	node_server();
	~node_server();
	node_server(const node_location&, const node_client& parent_id, real);
	integer get_boundary_size(std::vector<std::array<integer, NDIM>>&, std::vector<std::array<integer, NDIM>>&, integer,
			integer) const;
	void set_boundary(const std::vector<real>&, integer face);
	std::vector<real> get_boundary(integer face);

	integer regrid_gather();
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid_gather, regrid_gather_action);

	void regrid_scatter(integer, integer);
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid_scatter, regrid_scatter_action);

	void recv_boundary(const std::vector<real>&, integer rk, integer face);
	HPX_DEFINE_COMPONENT_ACTION(node_server, recv_boundary, send_boundary_action);

	real step();
	HPX_DEFINE_COMPONENT_ACTION(node_server, step, step_action);

	void regrid();
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid, regrid_action);

};

#endif /* NODE_SERVER_HPP_ */
