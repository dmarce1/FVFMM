#ifndef NDEBUG
#include <fenv.h>
#else
#include <mpi.h>
#endif

#include "defs.hpp"

#include "node_server.hpp"
#include "node_client.hpp"


HPX_PLAIN_ACTION(node_server::output, output_action_type);

void node_server::start_run() {

	printf("Starting...\n");
	solve_gravity(false, 3);

	double output_dt = 0.01;
	integer output_cnt = 0;

	real t = ZERO;
	integer step_num = 0;

	while (true) {

		if (t / output_dt >= output_cnt) {
			char* fname;
			if (asprintf(&fname, "X.%i.silo", int(output_cnt))) {
			}
			++output_cnt;
			hpx::async < output_action_type > (hpx::find_here(), fname).get();
			free(fname);
			printf("--- output ---\n");
		}
		double tstart = MPI_Wtime();
		real dt = step();
		printf("%i %e %e %e\n", int(step_num), double(t), double(dt), MPI_Wtime() - tstart);
		t += dt;
		++step_num;
	}
}

int hpx_main(int argc, char* argv[]) {
#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#endif
	node_client root_id = node_client::create(hpx::find_here());
//	auto root_id = hpx::new_ < node_server > (hpx::find_here()).get();
	node_client root_client(root_id);
	//root_client.register_(node_location()).get();
	for (integer l = 0; l < MAX_LEVEL; ++l) {
		root_client.regrid().get();
		printf("++++++++++++\n");
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

