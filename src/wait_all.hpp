/*
 * wait_all.hpp
 *
 *  Created on: Sep 1, 2015
 *      Author: dmarce1
 */



#ifndef WAIT_ALL_HPP_
#define WAIT_ALL_HPP_


#include <boost/chrono.hpp>

namespace fvfmm {

template<class Iterator1, class Iterator2>
void wait_all(Iterator1 ib, Iterator2 ie) {
	auto i = ib;
	while(i != ie ) {
		i->get();
		++i;
	}
}
}


#endif /* WAIT_ALL_HPP_ */
