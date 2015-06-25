/*
 * node_client.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"
#include <boost/serialization/list.hpp>

node_client::node_client(const hpx::id_type& id) :
		base_type(id) {
}

node_client::node_client() :
		base_type() {
}

hpx::future<integer> node_client::regrid_gather() const {
	return hpx::async<typename node_server::regrid_gather_action>(get_gid());
}



hpx::future<void> node_client::regrid_scatter(integer a, integer b) const {
	return hpx::async<typename node_server::regrid_scatter_action>(get_gid(), a, b);
}

hpx::future<void> node_client::register_(const node_location& nloc) const {
	return hpx::register_id_with_basename("node_location", get_gid(), nloc.unique_id());
}

hpx::future<void> node_client::unregister(const node_location& nloc) const {
	return hpx::unregister_id_with_basename("node_location", nloc.unique_id());
}

hpx::future<void> node_client::send_boundary(const std::vector<real> data, integer face) const {
	return hpx::async<typename node_server::send_boundary_action>(get_gid(), data, face);
}


hpx::future<void> node_client::step() const {
	return hpx::async<typename node_server::step_action>(get_gid());
}


hpx::future<void> node_client::regrid() const {
	return hpx::async<typename node_server::regrid_action>(get_gid());
}

