#ifndef CHANNEL_HPP_
#define CHANNEL_HPP_

#include "defs.hpp"

template<class T>
class channel {
private:
#ifndef NDEBUG
	bool set;
#endif
	hpx::future<T> future;
	hpx::promise<T> promise;
public:
	channel() {
#ifndef NDEBUG
		set = false;
#endif
		future = promise.get_future();
	}
	~channel() = default;
	channel(const channel&) = delete;
	channel(channel&& other ) {
		future = std::move(other.future);
		promise = std::move(other.promise);
	}
	channel& operator=(channel&& other ) {
		future = std::move(other.future);
		promise = std::move(other.promise);
		return *this;
	}

	template<class U>
	void set_value( U value ) {
		assert(!set);
		promise.set_value(value);
#ifndef NDEBUG
		set = true;
#endif
	}

	T get() {
		T data = future.get();
		promise.reset();
		future = promise.get_future();
#ifndef NDEBUG
		set = false;
#endif
		return data;
	}

};


#endif /* CHANNEL_HPP_ */
