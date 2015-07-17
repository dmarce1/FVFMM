/*
 * node_client.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"
#include <boost/serialization/list.hpp>

hpx::future<void> node_client::form_tree(const hpx::id_type& id1, const hpx::id_type& id2, const  std::vector<hpx::shared_future<hpx::id_type>>& ids) {
	return hpx::async<typename node_server::form_tree_action>(get_gid(),id1, id2, std::move(ids));
}
hpx::future<hpx::id_type> node_client::copy_to_locality(const hpx::id_type& id) const {
	return hpx::async<typename node_server::copy_to_locality_action>(get_gid(), id);
}


hpx::future<void> node_client::save( const std::string& str) const {
	return hpx::async<typename node_server::save_action>(get_gid(), str);
}

hpx::future<void> node_client::load( const std::string& str) {
	return hpx::async<typename node_server::load_action>(get_gid(), str);
}

node_client::node_client(const hpx::id_type& id) :
base_type(id) {
}


node_client& node_client::operator=(const hpx::id_type& id) {
	static_cast<base_type*>(this)->operator=(id);
	return *this;
}

node_client::node_client(hpx::future<hpx::id_type>&& id) :
base_type(std::move(id)) {
}


node_client& node_client::operator=(hpx::future<hpx::id_type>&& id) {
	static_cast<base_type*>(this)->operator=(std::move(id));
	return *this;
}

node_client::node_client(const hpx::shared_future<hpx::id_type>& id) :
base_type(id) {
}


node_client& node_client::operator=(const hpx::shared_future<hpx::id_type>& id) {
	static_cast<base_type*>(this)->operator=(id);
	return *this;
}

node_client::node_client() :
		base_type() {
}

hpx::future<void> node_client::solve_gravity(bool ene, integer c) const {
	return hpx::async<typename node_server::solve_gravity_action>(get_gid(), ene, c);
}

hpx::future<integer> node_client::regrid_gather() const {
	return hpx::async<typename node_server::regrid_gather_action>(get_gid());
}

hpx::future<void> node_client::regrid_scatter(integer a, integer b) const {
	return hpx::async<typename node_server::regrid_scatter_action>(get_gid(), a, b);
}
/*
hpx::future<void> node_client::register_(const node_location& nloc) const {
	return hpx::async([=]() {
		auto fut = hpx::register_id_with_basename("node_location", get_gid(), nloc.unique_id());
		if( !fut.get() ) {
			printf( "Failed to register node at location %s\n", nloc.to_str().c_str());
		}
	});
}*/

hpx::future<hpx::id_type> node_client::get_child_client(integer ci) {
	if (get_gid() != hpx::invalid_id) {
		return hpx::async<typename node_server::get_child_client_action>(get_gid(), ci);
	} else {
		return hpx::make_ready_future(hpx::invalid_id);
	}
}
/*
hpx::future<void> node_client::unregister(const node_location& nloc) const {
	return hpx::unregister_id_with_basename("node_location", nloc.unique_id());
}*/

hpx::future<void> node_client::send_hydro_boundary(const std::vector<real> data, integer rk, integer face) const {
	return hpx::async<typename node_server::send_hydro_boundary_action>(get_gid(), data, rk, face);
}

hpx::future<void> node_client::send_gravity_boundary(const std::vector<real> data, integer face, integer c) const {
	return hpx::async<typename node_server::send_gravity_boundary_action>(get_gid(), data, face, c);
}

hpx::future<void> node_client::send_gravity_multipoles(const multipole_pass_type& data, integer ci, integer c) const {
	return hpx::async<typename node_server::send_gravity_multipoles_action>(get_gid(), data, ci, c);
}

hpx::future<void> node_client::send_gravity_expansions(const expansion_pass_type& data, integer c) const {
	return hpx::async<typename node_server::send_gravity_expansions_action>(get_gid(), data, c);
}

hpx::future<real> node_client::step() const {
	return hpx::async<typename node_server::step_action>(get_gid());
}

hpx::future<void> node_client::start_run() const {
	return hpx::async<typename node_server::start_run_action>(get_gid());
}

hpx::future<void> node_client::regrid() const {
	return hpx::async<typename node_server::regrid_action>(get_gid());
}

