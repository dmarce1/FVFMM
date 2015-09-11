
#ifndef CHANNEL_HPP_
#define CHANNEL_HPP_

#include "defs.hpp"

#include <atomic>

template<class T>
class channel {
private:
	hpx::lcos::local::counting_semaphore signal;
	T data;
	std::atomic<bool> full;
public:
	channel() : full(false){
	}
	~channel() = default;
	channel(const channel&) = delete;
	channel(channel&& other ) = delete;
	channel& operator=(channel&& other ) = delete;

	template<class U>
	void set_value( U value ) {
		assert(!full);
		data = std::move(value);
		full = true;
		signal.signal();
	}

	T get() {
		signal.wait();
		assert(full);
		full = false;
		return data;
	}

};


#endif /* CHANNEL_HPP_ */

