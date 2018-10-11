#include <iostream>
#include <fstream>
#include "md_timer.h"

using namespace std;

void start_timer(struct MD_Timer * time) {
	gettimeofday(&time->clock_holder, NULL);
}
