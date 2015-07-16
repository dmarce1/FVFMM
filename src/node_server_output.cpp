/*
 * node_server_output.cpp
 *
 *  Created on: Jul 16, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"

HPX_PLAIN_ACTION(node_server::output_form, output_form_action);
HPX_PLAIN_ACTION(node_server::output_collect, output_collect_action);

grid::output_list_type node_server::olist;

void node_server::output_form() {

	const auto localities = hpx::find_all_localities();

	integer me_loc, child_r, child_l;

	me_loc = hpx::get_locality_id();

	child_r = (((me_loc + 1) << 1) + 1) - 1;
	child_l = (((me_loc + 1) << 1) + 0) - 1;
	child_r = child_r >= integer(localities.size()) ? -1 : child_r;
	child_l = child_l >= integer(localities.size()) ? -1 : child_l;

	std::list<hpx::future<void>> futs;

	if (child_r != -1) {
		futs.push_back(hpx::async < output_form_action > (localities[child_r]));
	}
	if (child_l != -1) {
		futs.push_back(hpx::async < output_form_action > (localities[child_l]));
	}

	{
		boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		olist = grid::output_list_type();
		for (auto i = local_node_list.begin(); i != local_node_list.end(); ++i) {
			if (!(*i)->is_refined) {
				grid::output_list_type this_list = (*i)->grid_ptr->get_output_list();
				grid::merge_output_lists(olist, std::move(this_list));
			}
		}
	}

	hpx::wait_all(futs.begin(), futs.end());
}

grid::output_list_type node_server::output_collect(const std::string& filename) {

	const auto localities = hpx::find_all_localities();

	integer me_loc, child_r, child_l;

	me_loc = hpx::get_locality_id();

	child_r = (((me_loc + 1) << 1) + 1) - 1;
	child_l = (((me_loc + 1) << 1) + 0) - 1;
	child_r = child_r >= integer(localities.size()) ? -1 : child_r;
	child_l = child_l >= integer(localities.size()) ? -1 : child_l;

	std::list < hpx::future < grid::output_list_type >> futs;

	if (child_r != -1) {
		futs.push_back(hpx::async < output_collect_action > (localities[child_r], filename));
	}
	if (child_l != -1) {
		futs.push_back(hpx::async < output_collect_action > (localities[child_l], filename));
	}

	{
		boost::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		for (auto i = futs.begin(); i != futs.end(); ++i) {
			grid::merge_output_lists(olist, i->get());
		}
	}

	if (hpx::get_locality_id() == 0) {
		grid::output(std::move(olist), filename.c_str());

	}
	return std::move(olist);
}
