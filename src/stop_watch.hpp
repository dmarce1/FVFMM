/*
 * stop_watch.hpp
 *
 *  Created on: Jun 9, 2015
 *      Author: dmarce1
 */

#ifndef STOP_WATCH_HPP_
#define STOP_WATCH_HPP_

class stop_watch {
private:
	clock_t t;
	clock_t t0;
	bool on;
	;
public:
	stop_watch() {
		t = clock_t(0);
		on = false;
	}
	void start() {
		if (!on) {
			t0 = clock();
			on = true;
		}

	}
	void stop() {
		if (on) {
			t += clock() - t0;
			on = false;
		}
	}
	double get() const {
		return double(t)/1000000.0;
	}
};

#endif /* STOP_WATCH_HPP_ */
