#ifndef MECHSYS_STOPWATCH_H
#define MECHSYS_STOPWATCH_H

// Std Lib
#include <stdio.h>        // for printf
#include <sys/time.h>     // for gettimeofday, getrusage
#include <sys/resource.h> // for getrusage
#include <time.h>         // for localtime

// MPI
#ifdef HAS_MPI
  #include <mpi.h>
#endif

// proc
#ifdef HAS_PROC
  #include <proc/readproc.h>
#endif

// mechsys
#include <mechsys/util/string.h>

#define CONVERT_TO_HMS(seconds,H,M,S)                         \
    H  = (unsigned int) ((unsigned int)seconds / 3600);       \
    M  = (unsigned int)(((unsigned int)seconds % 3600) / 60); \
    S  = seconds-(unsigned int)(H*3600)-(unsigned int)(M*60);

namespace Util
{

class Stopwatch
{
public:
    // Constructor
     Stopwatch (bool OnlyRoot=false, bool MemUsage=true);

    // Destructor
    ~Stopwatch (); ///< Will output time during destruction

    // Methods
    double CPUTime () const; ///< Returns the total CPU time consumed by the current process

private:
    bool    _only_root; ///< Only root prints timing (if MPI is activated)
    bool    _mem_usage; ///< Show memory usage as well ?
    timeval _start;     ///< Start time
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Stopwatch::Stopwatch (bool OnlyRoot, bool MemUsage)
    : _only_root(OnlyRoot), _mem_usage(MemUsage)
{
    // initial time
    gettimeofday (&_start, NULL);
}

inline Stopwatch::~Stopwatch ()
{
    // skip if not root
    if (_only_root)
    {
#ifdef HAS_MPI
        if (MPI::COMM_WORLD.Get_rank()!=0) return;
#endif
    }

    // final time
    timeval end;
    gettimeofday (&end, NULL);

    // interval
    double dif = end.tv_sec + (end.tv_usec/1000000.0) - _start.tv_sec - (_start.tv_usec/1000000.0);
    double cpu = CPUTime();

    // time H:M:S
    tm   t1 = (*(localtime (&_start.tv_sec)));
    tm * t2 =    localtime (&   end.tv_sec);

    // output time
    if (dif>60)
    {
        int h,m,H,M;
        double s,S;
        CONVERT_TO_HMS (dif,h,m,s)
        CONVERT_TO_HMS (cpu,H,M,S)
        printf("%s  Elapsed time       = %dh%dm%gs  CPU %dh%dm%gs  (%d:%02d:%02d => %d:%02d:%02d)%s\n", TERM_CLR3, h,m,s, H,M,S,
                t1. tm_hour, t1. tm_min, t1. tm_sec,
                t2->tm_hour, t2->tm_min, t2->tm_sec, TERM_RST);
    }
    else
    {
        printf("%s  Elapsed time       = %.6lfs  CPU %.6lfs  (%d:%02d:%02d => %d:%02d:%02d)%s\n", TERM_CLR3, dif, cpu,
                t1. tm_hour, t1. tm_min, t1. tm_sec,
                t2->tm_hour, t2->tm_min, t2->tm_sec, TERM_RST);
    }

    // memory usage
    if (_mem_usage)
    {
#ifdef HAS_PROC
        proc_t p;
        look_up_our_self (&p);
        printf("%s  Process memory     = %lu [kb]  %lu [Mb]%s\n", TERM_CLR5, p.vsize/1024, p.vsize/1048576, TERM_RST);
#endif
    }
}

inline double Stopwatch::CPUTime() const
{
    /* The timeval struct used to measure time has only two fields, and
       both are unsigned ints. They are named tv_sec and tv_usec, and
       jointly represent one single value. tv_sec*1000000+tv_usec gives the
       number of microseconds. http://rabbit.eng.miami.edu/info/functions/time.html#gtod */

    /* This function returns the total CPU time consumed by the current
       process, measured in seconds, as a double precision floating point number.
       It adds together the user time and the system time.
       Note: Although the format used is capable of measuring time to an accuracy
       of a microsecond, do not expect that much precision from any real system. */

    timeval tim;
    rusage  ru;
    getrusage (RUSAGE_SELF, &ru);
    tim = ru.ru_utime;
    double t = (double)tim.tv_sec + (double)tim.tv_usec/1000000.0;
    tim = ru.ru_stime;
    t += (double)tim.tv_sec + (double)tim.tv_usec/1000000.0;
    return t;
}

}; // namespace Util

#endif // MECHSYS_STOPWATCH_H
