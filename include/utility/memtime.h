/*
 * utility_general.h
 *
 *  Created on: Jun 17, 2018
 *      Author: ivan
 */

#ifndef SRC_UTILITY_MEMTIME_H_
#define SRC_UTILITY_MEMTIME_H_

#include <time.h>
#include <string>
#include <stdarg.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>

namespace raptor {

std::string GetUTCTime(std::string fmt="%a, %d %b %y %T %z");
std::string GetLocalTime();
// void ProcessMemUsage(double& vm_usage, double& resident_set);
std::string FormatMemoryConsumptionAsString();
size_t getCurrentRSS();
size_t getPeakRSS();

}

#endif
