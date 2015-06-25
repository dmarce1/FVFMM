/*
 * channel.hpp
 *
 *  Created on: Jun 19, 2015
 *      Author: dmarce1
 */

#ifndef CHANNEL_HPP_
#define CHANNEL_HPP_

#include "defs.hpp"

template<class T>
class channel {
private:
	bool set;
	hpx::future<T> future;
	hpx::promise<T> promise;
public:
	channel() {
		future = promise.get_future();
	}
	~channel() = default;
	channel(channel&&) = default;
	channel(const channel&) = delete;
	channel& operator=(channel&&) = default;
	channel& operator=(const channel&) = delete;

	template<class U>
	void set_value( U value ) {
		promise.set_value(value);
	}

	T get() {
		T data = future.get();
		promise.reset();
		future = promise.get_future();
		return data;
	}

};

template<>
class channel<void> {
private:
	bool set;
	hpx::shared_future<void> future;
	hpx::promise<void> promise;
public:
	channel();
	~channel() = default;
	channel(channel&&) = default;
	channel(const channel&) = delete;
	channel& operator=(channel&&) = default;
	channel& operator=(const channel&) = delete;

	void set_value();
	void get();
};

#endif /* CHANNEL_HPP_ */