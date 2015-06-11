#include "grid.hpp"
#include "problem.hpp"
#include "taylor.hpp"
#ifndef NDEBUG
#include <fenv.h>
#endif

int main(int argc, char* argv[]) {
	taylor<3> test;
#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#endif
	char* fname;
	grid root(star);
	real t = 0.0;
	integer step = 0;
	if (argc > 1) {
		root.load(argv[1]);
		t = root.get_time();
		step = root.get_step();
	}
	real tmax = 1000.0;
	while (t < tmax) {
		auto tstart = clock();
		real dt = root.step();
		auto tstop = clock();
		printf("%i %e %e %e\n", int(step), double(t), double(dt), double(tstop-tstart));
		t += dt;
		if (step % 10 == 0) {
			if (!asprintf(&fname, "X.%i.silo", int(step/10))) {
				abort();
			}
			root.output(fname);
			free(fname);
		}
		++step;
	}
}

