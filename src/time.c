#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


double GWHL_getTime(){
    struct timespec time;

    if( clock_gettime(CLOCK_MONOTONIC, &time) == -1 ) {
      perror( "clock gettime" );
      return NAN;
    }

    return (time.tv_sec) + (time.tv_nsec ) / 1.0e9;
}


