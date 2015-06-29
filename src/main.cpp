#include "node_server.hpp"
#include "node_client.hpp"
#ifndef NDEBUG
#include <fenv.h>
#else
#include <mpi.h>
#endif

HPX_PLAIN_ACTION(node_server::output, output_action_type);

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
	hpx::async < output_action_type > (hpx::find_here(), "X.0.silo").get();
	real t = ZERO;
	real tmax = 0.05;
	integer step_num = 0;
//	while (t < tmax) {
		real dt = root_client.step().get();
		printf("%i %e %e\n", int(step_num), double(t), double(dt));
		t += dt;
		++step_num;
//	}
	hpx::async < output_action_type > (hpx::find_here(), "X.1.silo").get();
	root_client.unregister(node_location()).get();
	printf("Exiting...\n");
	return hpx::finalize();
}

