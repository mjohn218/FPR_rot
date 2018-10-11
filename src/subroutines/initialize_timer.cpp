#include <iostream>
#include <fstream>
#include "md_timer.h"

using namespace std;

void initialize_timer(struct MD_Timer * time) {
	time->clock_holder.tv_sec = 0;
	time->clock_holder.tv_usec = 0;
	time->duration.tv_sec = 0;
	time->duration.tv_usec = 0;
}
