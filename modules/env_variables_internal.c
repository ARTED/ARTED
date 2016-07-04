/*
  Ab-initio Real-Time Electron Dynamics Simulator, ARTED
  Copyright (C) 2016  ARTED developers

  This file is part of env_variables_internal.c.

  env_variables_internal.c is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  env_variables_internal.c is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with env_variables_internal.c.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdlib.h>
#include <string.h>

#define ARTED_CPU_TASK_ENV      "ARTED_CPU_TASK_RATIO"       /* 0.1 ~ 1.0 */
#define ARTED_CPU_PPN_ENV       "ARTED_CPU_PPN"              /* Process per Node */
#define ARTED_MIC_PPN_ENV       "ARTED_MIC_PPN"              /* Process per Node */
#define ARTED_LOAD_BALANCER_ENV "ARTED_ENABLE_LOAD_BALANCER" /* 1 or 0 */

#define ARTED_SYM_DEBUG_ENV     "ARTED_SYM_DEBUG"            /* 1 or 0 */

void get_cpu_task_ratio_internal_(double * ret) {
  char* env = getenv(ARTED_CPU_TASK_ENV);
  if(env != NULL)
    *ret = atof(env);
  else
    *ret = 1.0;
}

void get_cpu_ppn_internal_(int * ret) {
  char* env = getenv(ARTED_CPU_PPN_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 1;
}

void get_mic_ppn_internal_(int * ret) {
  char* env = getenv(ARTED_MIC_PPN_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 1;
}

void get_load_balancer_flag_internal_(int * ret) {
  char* env = getenv(ARTED_LOAD_BALANCER_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 0;
}

void get_sym_debug_mode_internal_(int * ret) {
  char* env = getenv(ARTED_SYM_DEBUG_ENV);
  if(env != NULL)
    *ret = atoi(env);
  else
    *ret = 0;
}
