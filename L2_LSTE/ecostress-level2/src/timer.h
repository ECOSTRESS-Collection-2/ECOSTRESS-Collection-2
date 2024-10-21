// timer.h
#pragma once
#include <sys/time.h>
#ifndef NULL
#define NULL ((void*)0)
#endif
static time_t getusecs()
{
    struct timeval tmv;
    gettimeofday(&tmv, NULL);
    time_t usecs = tmv.tv_sec * 1000000 + tmv.tv_usec;
    return usecs;
}

static time_t elapsed(time_t since_time) {return getusecs() - since_time;}

static double toseconds(time_t usecs) {return (double)usecs * 0.000001;}
