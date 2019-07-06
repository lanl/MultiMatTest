/*
 * Copyright (c) 2017-2019, Triad National Security, LLC.
 * All rights Reserved.
 * 
 * This is the code released under LANL Copyright Disclosure C17041/LA-CC-17-041
 * Copyright 2017-2019.  Triad National Security, LLC. This material was produced
 * under U.S. Government contract 89233218CNA000001 for Los Alamos National
 * Laboratory (LANL), which is operated by Triad National Security, LLC
 * for the U.S. Department of Energy. See LICENSE file for details.
 *
 * Released under the New BSD License
 *
 * Bob Robey brobey@lanl.gov and Rao Garimella rao@lanl.gov
 */

#include <time.h>
#include "timer.h"

void cpu_timer_start(struct timespec *tstart_cpu)
{
   clock_gettime(CLOCK_MONOTONIC, tstart_cpu);
}
double cpu_timer_stop(struct timespec tstart_cpu)
{
   struct timespec tstop_cpu, tresult;
   clock_gettime(CLOCK_MONOTONIC, &tstop_cpu);
   tresult.tv_sec = tstop_cpu.tv_sec - tstart_cpu.tv_sec;
   tresult.tv_nsec = tstop_cpu.tv_nsec - tstart_cpu.tv_nsec;
   double result = (double)tresult.tv_sec + (double)tresult.tv_nsec*1.0e-9;

   return(result);
}

