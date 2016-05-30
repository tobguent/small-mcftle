#ifndef _TIMER_INCLUDE_ONCE_
#define _TIMER_INCLUDE_ONCE_

#include <iostream>
#include <iomanip>
#include <chrono>

class Timer
{
public:
	void tic() {
		t_start = std::chrono::high_resolution_clock::now();
	}
	void toc() {
		auto t_end = std::chrono::high_resolution_clock::now();

		std::cout << std::fixed << std::setprecision(2) 
			<< "Wall clock time passed: "
			<< std::chrono::duration<double, std::milli>(t_end - t_start).count()
			<< " ms\n";
	}

private:
	std::chrono::high_resolution_clock::time_point t_start;
};

#endif