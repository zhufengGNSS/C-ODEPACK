/*
  Copyright (C) 2011  Akshay Srinivasan 
  <akshay@ncbs.res.in>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "c-odepack.h"

void simple_pendulum_field(double *qdot,
			   const double t, const double *q, void *data){
  double *parms;
  parms = (double *) data;

  double alpha;
  alpha = parms[0];

  qdot[0] = q[1];
  qdot[1] = - alpha * sin(q[0]);

  return;
}

void simple_pendulum_jacobian(double *dfdq,
			      const double t, const double *q, void *data){
  double *parms;
  parms = (double *) data;

  double alpha;
  alpha = parms[0];
  
/*Column ordered; ODEPACK is in FORTRAN 77.*/
  dfdq[0] = 0.0; dfdq[2] = 1.0;
  
  dfdq[1] = - alpha * cos(q[0]); dfdq[3] = 0.0;

  return;
}

void simple_pendulum_cycle(double *g,
			   const double t, const double *q, void *data){
  double *parms;
  parms = (double *) data;
  
  double alpha;
  alpha = parms[0];

  double theta_0;
  theta_0 = parms[1];
/*-----------------------------------------------------------------------*/
  g[0] = (q[0] - theta_0) * q[1];

  /*Make sure the initial value is not detected!*/
  if(t < 1e-2)
    g[0] += 0.1;

  return;
}

int main(){
/*-----------------------------------------------------------------------*/  
  double opkd_rtol = 0.0, opkd_atol = 1e-12;
  dlsodar_options *opkd_opt;
  opkd_opt = calloc(1, sizeof(dlsodar_options));
  
  opkd_opt->NEQ = 2;
  opkd_opt->NG = 1;
  
  opkd_opt->jac_type = USR_FULL_JAC;
  
  opkd_opt->max_steps = 10000;
  opkd_opt->atol = &opkd_atol;
  opkd_opt->rtol = &opkd_rtol;

  dlsodar_workspace *opkd_work;
  opkd_work = dlsodar_workspace_alloc(opkd_opt);
/*-----------------------------------------------------------------------*/
  double q[2], t0, tf, dt, t;

  t0 = 0.0;
  tf = 100.0;
  dt = 0.1;

  q[0] = M_PI * 999/1000.0;
  q[1] = 0.;

  double alpha = 1.0;
  double parms[2] = {alpha, q[0]};
/*-----------------------------------------------------------------------*/  
  FILE *data;
  data = fopen("simple_pendulum_trajectory", "w");

  int i = 0, cyc = 0;
  t = t0;
  
  while(t < tf){    
    dlsodar_integrate(t + dt,
		      &t, q,
		      &simple_pendulum_field, NULL// &simple_pendulum_jacobian
		      , &simple_pendulum_cycle,
		      &parms, opkd_work);
    fprintf(data, "%.15lf\t%.15lf\t%.15lf\n", t + dt, q[0], q[1]);
    
    if(opkd_work->jroot[0] == 1)
      cyc += 1;
    if(cyc == 2)
      break;
  }

  fclose(data);
/*-----------------------------------------------------------------------*/  
  dlsodar_workspace_free(opkd_work);
  free(opkd_opt);
  return;
}
