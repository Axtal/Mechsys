#ifndef MECHSYS_STOPWATCH_H
#define MECHSYS_STOPWATCH_H

// Std Lib
#include <stdio.h>        // for printf
#include <sys/time.h>     // for gettimeofday, getrusage
#include <sys/resource.h> // for getrusage
#include <time.h>         // for localtime

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
    Stopwatch (bool MemUsage=true) : _mem(MemUsage)
    {
        gettimeofday (&_start, NULL);
    }

    ~Stopwatch ()
    {
        timeval end;
        gettimeofday (&end, NULL);

        /*
        long seconds  = end.tv_sec  - _start.tv_sec;
        long useconds = end.tv_usec - _start.tv_usec;
        long mtime    = ((seconds) * 1000 + useconds/1000.0) + 0.5;
        printf("milliseconds = %ld ms",mtime);
        */

        double d_sta = _start.tv_sec + (_start.tv_usec/1000000.0);
        double d_end =    end.tv_sec +    (end.tv_usec/1000000.0);

        tm * t1 = localtime (&_start.tv_sec);
        tm * t2 = localtime (&   end.tv_sec);
        printf("%s  Elapsed time       = %.6lf s  CPU %.6lf s  (%d:%02d:%02d => %d:%02d:%02d)%s\n", TERM_CLR3,
                d_end-d_sta, CPUTime(),
                t1->tm_hour, t1->tm_min, t1->tm_sec,
                t2->tm_hour, t2->tm_min, t2->tm_sec, TERM_RST);

        if (_mem)
        {
#ifdef HAS_PROC
            proc_t p;
            look_up_our_self (&p);
            printf("%s  Process memory     = %lu [kb]  %lu [Mb]%s\n", TERM_CLR5, p.vsize/1024, p.vsize/1048576, TERM_RST);
#endif
        }
    }

    /* The timeval struct used to measure time has only two fields, and
       both are unsigned ints. They are named tv_sec and tv_usec, and
       jointly represent one single value. tv_sec*1000000+tv_usec gives the
       number of microseconds. http://rabbit.eng.miami.edu/info/functions/time.html#gtod */

    double CPUTime ()
    {
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

private:
    bool    _mem;
    timeval _start;
    double  _d_start;
};

}; // namespace Util

#endif // MECHSYS_STOPWATCH_H
