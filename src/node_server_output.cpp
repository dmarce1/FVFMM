/*
 * node_server_output.cpp
 *
 *  Created on: Jul 16, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"

HPX_PLAIN_ACTION(node_server::output_form, output_form_action);
HPX_PLAIN_ACTION(node_server::output_collect, output_collect_action);

std::list<grid::output_list_type> node_server::olists;

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
		for (auto i = local_node_list.begin(); i != local_node_list.end(); ++i) {
			if (!(*i)->is_refined) {
				olists.push_back((*i)->grid_ptr->get_output_list());
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

	grid::output_list_type olist;

	while (olists.size() > 1) {
	//	printf( "." );
		fflush(stdout);
		std::list<hpx::future<void>> merge_futs;
		std::vector<grid::output_list_type> v(olists.size());
		auto j = olists.begin();
		for (integer i = 0; i != integer(olists.size()); ++i) {
			v[i] = std::move(*j);
			++j;
		}
		olists.clear();
		for (integer i = 0; i < integer(v.size()); i += 2) {
			if (i == integer(v.size()) - 1) {
				olists.push_back(std::move(v[i]));
			} else {
				merge_futs.push_back( hpx::async([&](integer i1, integer i2)-> void {
					auto list1 = std::move(v[i1]);
					auto list2 = std::move(v[i2]);
					grid::merge_output_lists(list1, list2);
					olists.push_back(std::move(list1));
				}, i, i + 1));
			}
		}
		hpx::wait_all(merge_futs.begin(), merge_futs.end());
	}

	olist = std::move(olists.front());
	olists.clear();
	for (auto i = futs.begin(); i != futs.end(); ++i) {
//		printf( "+" );
		fflush(stdout);
		auto tmp = i->get();
		grid::merge_output_lists(olist, tmp);
	}
//	printf( "*" );
	fflush(stdout);

	if (hpx::get_locality_id() == 0) {
	//	printf( "1\n" );
		grid::output(olist, filename.c_str());
	//	printf( "0\n" );
	}
	return std::move(olist);
}


void node_server::save(std::string const&) const {

}

void node_server::load(std::string const&)  {

}
