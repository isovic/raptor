/*
 * utility_general.cc
 *
 *  Created on: Oct 11, 2014
 *      Author: ivan
 */

#include <unistd.h>
#include <ios>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <dirent.h>
#include <utility/memtime.h>

namespace raptor {

std::string GetUTCTime(std::string fmt) {
  //  const char* fmt = "%a, %d %b %y %T %z";
  char outstr[200];
  time_t t;
  struct tm *tmp;

  t = time(NULL);
  tmp = gmtime(&t);

  if (tmp == NULL) {
    fprintf (stderr, "ERROR: gmtime returned with error!\n");
    fflush (stderr);
    return std::string("[no_time]");
  }

  if (strftime(outstr, sizeof(outstr), fmt.c_str(), tmp) == 0) {
    fprintf (stderr, "ERROR: Problem formatting time into string!\n");
    fflush (stderr);
    return std::string("[no_time]");
  }

  return std::string(outstr);
}

std::string GetLocalTime() {
  char outstr[200];

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  sprintf(outstr, "%02d:%02d:%02d", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

  return std::string(outstr);
}



// //////////////////////////////////////////////////////////////////////////////
// //
// // process_mem_usage(double &, double &) - takes two doubles by reference,
// // attempts to read the system-dependent data for a process' virtual memory
// // size and resident set size, and return the results in KB.
// //
// // On failure, returns 0.0, 0.0
// //
// // Source: http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c

// void ProcessMemUsage(double& vm_usage, double& resident_set)
// {
//    using std::ios_base;
//    using std::ifstream;
//    using std::string;

//    vm_usage     = 0.0;
//    resident_set = 0.0;

//    // 'file' stat seems to give the most reliable results
//    //
//    ifstream stat_stream("/proc/self/stat",ios_base::in);

//    // dummy vars for leading entries in stat that we don't care about
//    //
//    string pid, comm, state, ppid, pgrp, session, tty_nr;
//    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
//    string utime, stime, cutime, cstime, priority, nice;
//    string O, itrealvalue, starttime;

//    // the two fields we want
//    //
//    unsigned long vsize;
//    long rss;

//    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
//                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
//                >> utime >> stime >> cutime >> cstime >> priority >> nice
//                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

//    stat_stream.close();

//    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
//    vm_usage     = vsize / 1024.0;
//    resident_set = rss * page_size_kb;
// }

// std::string FormatMemoryConsumptionAsString() {
//   std::stringstream ss;
//   ss << "[currentRSS = " << getCurrentRSS()/(1024*1024) << " MB, peakRSS = " << getPeakRSS()/(1024*1024) << " MB]";
//   return ss.str();
// }



/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}

}
