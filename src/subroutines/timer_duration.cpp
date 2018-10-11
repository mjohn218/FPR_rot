#include <iostream>
#include <fstream>
#include "md_timer.h"

using namespace std;

double timer_duration(const struct MD_Timer time) {
	return time.duration.tv_sec + 1.0e-6 * (double) time.duration.tv_usec;
}
