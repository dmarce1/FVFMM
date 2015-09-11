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
	std::list<std::shared_ptr<hpx::promise<T>>> pq;
	std::list<std::shared_ptr<hpx::future<T>>> fq;
	using mutex_type = 	hpx::lcos::local::spinlock;
	mutex_type mtx;

	void add_promise() {
		auto p_ptr = std::make_shared<hpx::promise<T>>();
		auto f_ptr = std::make_shared<hpx::future<T>>(p_ptr->get_future());
		pq.push_back(p_ptr);
		fq.push_back(f_ptr);
	}

public:
	channel() {
		boost::lock_guard<mutex_type> lock(mtx);
		add_promise();
	}
	~channel() = default;
	channel(const channel&) = delete;
	channel(channel&& other ) = delete;
	channel& operator=(channel&& other ) = delete;

	void set_value( T value ) {
		std::shared_ptr<hpx::promise<T>> p_ptr;
		{
			boost::lock_guard<mutex_type> lock(mtx);
			p_ptr = *(pq.begin());
			pq.pop_front();
		}
		p_ptr->set_value(std::move(value));
	}

	T get() {
		std::shared_ptr<hpx::future<T>> f_ptr;
		{
			boost::unique_lock<mutex_type> lock(mtx);
			f_ptr = *(fq.begin());
			fq.pop_front();
			add_promise();
		}
		return f_ptr->get();
	}

};


#endif /* CHANNEL_HPP_ */
