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
    gettimeofday (&_start, NULL);
}

inline Stopwatch::~Stopwatch ()
{
    if (_only_root)
    {
#ifdef HAS_MPI
        if (MPI::COMM_WORLD.Get_rank()!=0) return;
#endif
    }

    timeval end;
    gettimeofday (&end, NULL);

    double d_sta = _start.tv_sec + (_start.tv_usec/1000000.0);
    double d_end =    end.tv_sec +    (end.tv_usec/1000000.0);

    // This is the only way to get different values from local time. TODO: improve this
    tm * t = localtime (&_start.tv_sec);
    char t1[32];
    sprintf (t1, "%d:%02d:%02d", t->tm_hour, t->tm_min, t->tm_sec);

    t = localtime (&end.tv_sec);
    char t2[32];
    sprintf (t2, "%d:%02d:%02d", t->tm_hour, t->tm_min, t->tm_sec);

    printf("%s  Elapsed time       = %.6lf s  CPU %.6lf s  (%s => %s)%s\n", TERM_CLR3, d_end-d_sta, CPUTime(), t1, t2, TERM_RST);

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
