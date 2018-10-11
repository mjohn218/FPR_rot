#include <iostream>
#include <fstream>
#include "md_timer.h"

using namespace std;

void stop_timer(struct MD_Timer * time) {
	struct timeval end_tv;
	gettimeofday(&end_tv, NULL);
	time->duration.tv_sec += (end_tv.tv_sec - time->clock_holder.tv_sec);
	time->duration.tv_usec += (end_tv.tv_usec - time->clock_holder.tv_usec);
}
