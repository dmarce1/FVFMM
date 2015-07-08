#include "node_server.hpp"
#include "node_client.hpp"
#ifndef NDEBUG
#include <fenv.h>
#else
#include <mpi.h>
#endif

HPX_PLAIN_ACTION(node_server::output, output_action_type);

void node_server::start_run() {

	solve_gravity(false);

	double output_dt = 0.1;
	integer output_cnt = 0;

	real t = ZERO;
	integer step_num = 0;
	while (true) {

		if (t / output_dt >= output_cnt) {
			char* fname;
			if (asprintf(&fname, "X.%i.silo", int(output_cnt))) {
			}
			++output_cnt;
	//		hpx::async < output_action_type > (hpx::find_here(), fname).get();
			free(fname);
			printf("--- output ---\n");
		}
		real dt = step();
		printf("%i %e %e\n", int(step_num), double(t), double(dt));
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
	auto root_id = hpx::new_ < node_server > (hpx::find_here()).get();
	node_client root_client(root_id);
	root_client.register_(node_location()).get();
	for (integer l = 0; l < MAX_LEVEL; ++l) {
		root_client.regrid().get();
	}

	root_client.start_run().get();

	root_client.unregister(node_location()).get();
	printf("Exiting...\n");
	return hpx::finalize();
}

