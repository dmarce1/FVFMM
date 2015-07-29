/*
 * node_server_output.cpp
 *
 *  Created on: Jul 16, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"
#include <mutex>
#include <thread>

HPX_PLAIN_ACTION(node_server::output_form, output_form_action);
HPX_PLAIN_ACTION(node_server::output_collect, output_collect_action);

HPX_REGISTER_PLAIN_ACTION( save_action);

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
		std::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
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



std::size_t node_server::load_me( FILE* fp, integer& locality_id) {
	std::size_t cnt = 0;
	auto foo = std::fread;

//cnt += foo(&is_refined, sizeof(bool), 1, fp)*sizeof(bool);
	cnt += foo(&locality_id, sizeof(integer), 1, fp)*sizeof(integer);
	cnt += foo(&step_num, sizeof(integer), 1, fp)*sizeof(integer);
	cnt += foo(&current_time,sizeof(real), 1, fp)*sizeof(real);
	cnt += foo(&dx, sizeof(real), 1,  fp)*sizeof(real);
	cnt += foo(&(xmin[0]), sizeof(real), NDIM, fp)*sizeof(real);

	cnt += my_location.load(fp);
	cnt += grid_ptr->load(fp);
	return cnt;
}

std::size_t node_server::save_me( FILE* fp ) const {
	auto foo = std::fwrite;
	std::size_t cnt = 0;

//	cnt += foo(&is_refined, sizeof(bool), 1, fp)*sizeof(bool);
	integer locality_id = hpx::get_locality_id();
	cnt += foo(&locality_id, sizeof(integer), 1, fp)*sizeof(integer);
	cnt += foo(&step_num, sizeof(integer), 1, fp)*sizeof(integer);
	cnt += foo(&current_time,sizeof(real), 1, fp)*sizeof(real);
	cnt += foo(&dx, sizeof(real), 1, fp)*sizeof(real);
	cnt += foo(&(xmin[0]), sizeof(real), NDIM, fp)*sizeof(real);

	cnt += my_location.save(fp);
	cnt += grid_ptr->save(fp);
	return cnt;
}

void my_system( const std::string& command) {
//	printf( "Executing system command: %s\n", command.c_str());
	if( system( command.c_str()) != EXIT_SUCCESS ) {
		assert(false);
	}
}

std::pair<integer,std::size_t> node_server::save(const std::string& filename) {
	integer cnt = 0;
	std::size_t bytes_written = 0;
	std::list<hpx::future<std::pair<integer,std::size_t>>> save_futs;
	const auto my_id = hpx::get_locality_id();
	const auto localities = hpx::find_all_localities();
	if( my_id == 0 ) {
		for( auto i = localities.begin(); i != localities.end(); ++i) {
			if( *i != hpx::find_here() ) {
				save_futs.push_back(hpx::async<save_action>(*i, filename));
			}
		}
	}

	std::string this_name = filename + "." + std::to_string(integer(my_id));
	FILE* fp = fopen(this_name.c_str(), "wb");
	{
		std::lock_guard<hpx::lcos::local::spinlock> lock(local_node_list_lock);
		for (auto i = local_node_list.begin(); i != local_node_list.end(); ++i) {
//			printf( "Saving at %s\n",(*i)->my_location.to_str().c_str());
			bytes_written += (*i)->my_location.save(fp);
			bytes_written += (*i)->save_me(fp);
			++cnt;
		}
	}
	fclose(fp);
	std::size_t total_cnt = cnt;
	if( my_id == 0 ) {
		for( auto i = save_futs.begin(); i != save_futs.end(); ++i) {
			auto tmp =  i->get();
			total_cnt += tmp.first;
			bytes_written += tmp.second;
		}
		FILE* fp = fopen("size.tmp2", "wb");
		bytes_written += 2 * fwrite( &total_cnt, sizeof(integer), 1, fp) * sizeof(integer);
		std::size_t tmp = bytes_written + 2 * sizeof(std::size_t) + sizeof(real);
		bytes_written += 2 * fwrite( &tmp, sizeof(std::size_t), 1, fp)*sizeof(std::size_t);
		fclose( fp);
		my_system( "cp size.tmp2 size.tmp1\n");
		fp = fopen("size.tmp1", "ab");
		real omega = grid::get_omega();
		bytes_written += fwrite( &omega, sizeof(real), 1, fp) * sizeof(real);
		fclose( fp);
		std::string command = "rm -r -f " + filename + "\n";
		my_system (command);
		command = "cat size.tmp1 ";
			for( std::size_t i = 0; i != localities.size(); ++i ){
			command += filename + "." + std::to_string(integer(i)) + " ";
		}
		command += "size.tmp2 > " + filename + "\n";
		my_system( command);
		command = "rm -f -r " + filename + ".*\n";
		my_system( command);
		my_system( "rm -r -f size.tmp\n");
		printf( "Saved %i sub-grids with %lli bytes written\n", int(total_cnt), (long long)(bytes_written));

	}
	return std::make_pair(total_cnt, bytes_written);
}

void node_server::load(const std::string& filename, node_client root) {
	FILE* fp = fopen( filename.c_str(), "rb");
	integer cnt = 0;
	std::size_t bytes_expected;
	std::size_t bytes_read = 0;
	integer total_cnt;
	real omega;
	bytes_read += fread(&total_cnt, sizeof(integer), 1, fp)*sizeof(integer);
	bytes_read += fread(&bytes_expected, sizeof(std::size_t), 1, fp)*sizeof(std::size_t);
	bytes_read += fread(&omega, sizeof(real), 1, fp)*sizeof(real);
	grid::set_omega(omega);
	printf( "Loading %i subgrids...\n", int(total_cnt));
	for( integer i = 0; i != total_cnt; ++i ) {
		auto ns = std::make_shared<node_server>();
		node_location next_loc;
		bytes_read += next_loc.load(fp);
		std::size_t file_pos = bytes_read;
		integer dummy;
		bytes_read += ns->load_me(fp, dummy);
//		printf( "Loading at %s\n", next_loc.to_str().c_str());
		root.load_node(file_pos, filename,  next_loc, root.get_gid()).get();
		++cnt;
	}
	std::size_t bytes_check;
	integer ngrid_check;
	bytes_read += fread(&ngrid_check, sizeof(integer), 1, fp)*sizeof(integer);
	bytes_read += fread(&bytes_check, sizeof(std::size_t), 1, fp)*sizeof(std::size_t);
	fclose(fp);
	if( bytes_expected != bytes_read || bytes_expected != bytes_check || ngrid_check != cnt ) {
		printf( "Checkpoint file corrupt\n");
		printf( "Bytes read        = %lli\n", (long long) bytes_read);
		printf( "Bytes expected    = %lli\n", (long long) bytes_expected);
		printf( "Bytes end         = %lli\n", (long long) bytes_check);
		printf( "subgrids expected = %i\n", (int) total_cnt);
		printf( "subgrids end      = %i\n", (int) ngrid_check);
		abort();
	} else {
		printf( "--------Loaded: %i subgrids, %lli bytes read----------\n",int(cnt), (long long)(bytes_read));
	}
}



hpx::id_type node_server::load_node(std::size_t filepos, const std::string& fname, const node_location& loc, const hpx::id_type& _me) {
	me = _me;
	hpx::future<hpx::id_type> rc;
	if( loc == my_location) {
		integer locality_id;
		FILE* fp = fopen( fname.c_str(), "rb");
		fseek(fp, filepos, SEEK_SET);
		load_me(fp, locality_id);
		fclose(fp);
		if( locality_id != hpx::get_locality_id()) {
			const auto localities = hpx::find_all_localities();
			printf( "Moving %s from %i to %i\n", loc.to_str().c_str(), hpx::get_locality_id(), int(locality_id));
			rc =  me.copy_to_locality(localities[locality_id]);
		} else {
			rc = hpx::make_ready_future(me.get_gid());
		}
	} else {
		rc = hpx::make_ready_future(me.get_gid());
		if( !is_refined) {
			if( my_location.level() >= MAX_LEVEL) {
				abort();
			}
			children.resize(NCHILD);
			for (integer ci = 0; ci != NCHILD; ++ci) {
				children[ci] = hpx::new_<node_server>(hpx::find_here(), my_location.get_child(ci), me, ZERO);
			}
			is_refined = true;
			const integer flags = (my_location.level() == 0) ? GRID_IS_ROOT : 0;
			grid_ptr = std::make_shared < grid > (dx, xmin, flags);
		}
		for( integer ci = 0; ci != NCHILD; ++ci) {
			auto cloc =  my_location.get_child(ci);
			if( cloc == loc || loc.is_child_of(cloc)) {
				children[ci] = children[ci].load_node(filepos, fname, loc, children[ci].get_gid());
			}
		}
	}
	clear_family();
	return rc.get();

}

