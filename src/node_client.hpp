/*
 * node_client.hpp
 *
 *  Created on: Jun 11, 2015
 *      Author: dmarce1
 */

#ifndef NODE_CLIENT_HPP_
#define NODE_CLIENT_HPP_

#include "defs.hpp"
#include "node_location.hpp"

class node_server;

class node_client: public hpx::components::client_base<node_client, node_server> {
private:
	typedef hpx::components::client_base<node_client, node_server> base_type;
public:
	node_client();
	node_client(const hpx::id_type& id);
	hpx::future<void> regrid_scatter(integer, integer) const;
	hpx::future<void> register_(const node_location&) const;
	hpx::future<void> unregister(const node_location&) const;
	hpx::future<integer> regrid_gather() const;
	hpx::future<void> send_boundary(const std::vector<real>, integer rk, integer face) const;
	hpx::future<real> step() const;
	hpx::future<void> regrid() const;

};
#endif /* NODE_CLIENT_HPP_ */
