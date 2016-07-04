/*
  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
  Copyright (C) 2016  ARTED developers

  This file is part of papi_wrap.c.

  papi_wrap.c is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  papi_wrap.c is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with papi_wrap.c.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef ARTED_USE_PAPI

#include <papi.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

int    *EventSet;
double values[2];

void papi_begin_() {
  int ret, i;

  if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
    fprintf(stderr, "PAPI library init error!\n");
    exit(1);
  }
  if(PAPI_thread_init((unsigned long (*)(void))(omp_get_num_threads)) != PAPI_OK) {
    fprintf(stderr, "PAPI thread init error.\n");
    exit(1);
  }
  if (PAPI_num_counters() < 2) {
    fprintf(stderr, "No hardware counters here, or PAPI not supported.\n");
    exit(1);
  }

  EventSet = (int*) malloc(sizeof(int) * omp_get_max_threads());
  for(i = 0 ; i < omp_get_max_threads() ; EventSet[i++] = PAPI_NULL);

#pragma omp parallel
  {
    int t = omp_get_thread_num();
    PAPI_create_eventset(&EventSet[t]);
    PAPI_add_event(EventSet[t], PAPI_SP_OPS);
    PAPI_add_event(EventSet[t], PAPI_DP_OPS);

    if ((ret = PAPI_start(EventSet[t])) != PAPI_OK) {
      fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
      exit(1);
    }
  }

  for(i = 0 ; i < 2 ; values[i++] = 0);
}

void papi_end_() {
  int ret, i;
  long long v[2];
  long long v0, v1;
  double vin[2];

  v0 = v1 = 0;
#pragma omp parallel shared(v) reduction(+:v0,v1)
  {
    int t = omp_get_thread_num();
    if ((ret = PAPI_stop(EventSet[t],v)) != PAPI_OK) {
      fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
      exit(1);
    }
    v0 += v[0];
    v1 += v[1];
  }
  vin[0] = v0;
  vin[1] = v1;

  MPI_Reduce(vin, values, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  PAPI_shutdown();

  free(EventSet);
}

void papi_result_(double *time) {
  printf("SP FLOP = %f\n", values[0]);
  printf("DP FLOP = %f\n", values[1]);
  printf("GFLOPS  = %.2f\n", ((values[0] + values[1]) / *time) * 1.0e-9);
}

#else

void papi_begin_()              {}
void papi_end_()                {}
void papi_result_(double* time) {}

#endif

