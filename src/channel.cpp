/*
 * channel.cpp
 *
 *  Created on: Jun 19, 2015
 *      Author: dmarce1
 */

#include "channel.hpp"

channel<void>::channel() {
	future = promise.get_future();
}

void channel<void>::set_value() {
	promise.set_value();
	promise.reset();
	future = promise.get_future();
}

void channel<void>::get() {
	future.get();
}

