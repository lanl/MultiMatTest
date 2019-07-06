/*
 * Copyright (c) 2017, Los Alamos National Security, LLC.
 * All rights Reserved.
 * 
 * This is the code released under LANL Copyright Disclosure C17041/LA-CC-17-041
 * Copyright 2017.  Los Alamos National Security, LLC. This material was produced
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC
 * for the U.S. Department of Energy. See LICENSE file for details.
 *
 * Released under the New BSD License
 *
 * Bob Robey brobey@lanl.gov and Rao Garimella rao@lanl.gov
 */

#ifndef _TIMER_H
#define _TIMER_H

#ifdef __cplusplus
extern "C"
{
#endif

void cpu_timer_start(struct timespec *tstart_cpu);
double cpu_timer_stop(struct timespec tstart_cpu);

#ifdef __cplusplus
}
#endif

#endif  /* _TIMER_H */

