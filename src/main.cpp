#include <fenv.h>
#ifdef NDEBUG
#include <mpi.h>
#endif

#include "defs.hpp"

#include "node_server.hpp"
#include "node_client.hpp"

HPX_PLAIN_ACTION(node_server::output_collect, output_collect_action_type);
HPX_PLAIN_ACTION(node_server::output_form, output_form_action_type);

void node_server::start_run() {

	printf("Starting...\n");
	solve_gravity(false, 3);

	double output_dt = (real(2)*real(M_PI)/DEFAULT_OMEGA)/256.0;
	integer output_cnt;

	real& t = current_time;
	integer step_num = 0;

	node_server* root_ptr = me.get_ptr().get();

	hpx::future<void> output_done = hpx::make_ready_future();
	output_cnt = root_ptr->get_time() / output_dt;
	while (true) {

		auto diags = diagnostics();
		FILE* fp = fopen( "diag.dat", "at");
		fprintf( fp, "%12.6e ", t);
		for( integer f = 0; f != NF; ++f) {
			fprintf( fp, "%15.8e ", diags.grid_sum[f]+diags.outflow_sum[f]);
			fprintf( fp, "%15.8e ", diags.outflow_sum[f]);
			}
		fprintf( fp, "\n");
		fclose(fp);

		fp = fopen( "spin.dat", "at");
		fprintf( fp, "%12.6e ", t);
		for( integer f = 0; f != NDIM; ++f) {
			fprintf( fp, "%15.8e ", diags.s_sum[f]);
			fprintf( fp, "%15.8e ", diags.l_sum[f]);
			}
		fprintf( fp, "\n");
		fclose(fp);

#ifdef NDEBUG
		double tstart = MPI_Wtime();
#endif
		if (t / output_dt >= output_cnt) {
			output_done.get();
			char* fname;

			printf("--- begin checkpoint ---\n");
			if (asprintf(&fname, "X.%i.chk", int(output_cnt))) {
			}
			save_action()(hpx::find_all_localities()[0], fname);
			printf("--- end checkpoint ---\n");
			free(fname);

			if (asprintf(&fname, "X.%i.silo", int(output_cnt))) {
			}
			++output_cnt;
			printf("--- begin output ---\n");
			hpx::async < output_form_action_type > (hpx::find_here()).get();
			printf("--- middle output ---\n");
			std::string _fname = fname;
			output_done = (hpx::async < output_collect_action_type > (hpx::find_here(), std::move(_fname))).then(
					[](hpx::future<grid::output_list_type>&& fut) -> void {
						fut.get();
						printf("--- end output ---\n");
					});

			free(fname);

			if (asprintf(&fname, "X.%i.silo", int(output_cnt))) {
			}
			free(fname);
		}
		real dt = step();
#ifdef NDEBUG
		printf("%i %e %e %e\n", int(step_num), double(t), double(dt), MPI_Wtime() - tstart);
#else
		printf("%i %e %e\n", int(step_num), double(t), double(dt));
#endif
		t += dt;
		++step_num;

	}
}

int hpx_main(int argc, char* argv[]) {
//#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
//#endif
	node_client root_id = node_client::create(hpx::find_here());
	node_client root_client(root_id);

	if (argc == 1) {
		for (integer l = 0; l < MAX_LEVEL; ++l) {
			root_client.regrid().get();
			printf("++++++++++++\n");
		}
	} else {
		std::string fname(argv[1]);
		printf( "Loading from %s...\n", fname.c_str());
		root_client.get_ptr().get()->load(fname, root_client);
		root_client.regrid().get();
		printf( "Done. \n");
	}

	std::vector < hpx::shared_future < hpx::id_type >> null_sibs(NFACE);
	for (integer si = 0; si != NFACE; ++si) {
		null_sibs[si] = hpx::make_ready_future(hpx::invalid_id).share();
	}
	root_client.form_tree(root_client.get_gid(), hpx::invalid_id, null_sibs).get();
	printf("------------\n");

	root_client.start_run().get();

	//root_client.unregister(node_location()).get();
	printf("Exiting...\n");
	return hpx::finalize();
}

