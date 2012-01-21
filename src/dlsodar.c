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

#include "opkdmain.h"
#include "c-odepack.h"
#include "utility.h"

dlsodar_workspace * dlsodar_workspace_alloc(const dlsodar_options *opt){
  dlsodar_workspace *dls;
  dls = (dlsodar_workspace *) calloc(1, sizeof(dlsodar_workspace));
  did_malloc_work(dls, "dlsodar_workspace_alloc");
  
  dls->neq = opt->NEQ;
  dls->ng = opt->NG;
/*-----------------------------------------------------------------------*/
  dls->jroot = calloc(dls->ng, sizeof(int));
  did_malloc_work(dls->jroot, "dlsodar_workspace_alloc");
/*--------------------------------------------------------------*/  
  dls->lrw = 22 + dls->neq * max(16, dls->neq + 9) + 3 * dls->ng;
  dls->rwork = calloc(dls->lrw, sizeof(double));
  did_malloc_work(dls->rwork, "dlsodar_workspace_alloc");		  
/*--------------------------------------------------------------*/
  dls->liw = 20 + dls->neq;
  dls->iwork = calloc(dls->liw, sizeof(int));
  did_malloc_work(dls->iwork, "dlsodar_workspace_alloc");
/*-----------------------------------------------------------------------*/  
  dls->iopt = 1;     

  dls->rwork[4] = opt->init_step_size;
  dls->rwork[5] = opt->max_step_size;
  dls->rwork[6] = opt->min_step_size;

  dls->iwork[4] = opt->print_at_switch;
  dls->iwork[5] = opt->max_steps;
  dls->iwork[6] = opt->max_num_messages;
  dls->iwork[7] = opt->max_order_non_stiff;
  dls->iwork[8] = opt->max_order_stiff;
/*--------------------------------------------------------------*/
  if((opt->itol < 2) || (opt->itol > 4)){
    dls->itol = 1;//Error bounds are same for all dimensions.
  }
  else
    dls->itol = opt->itol;
  
  dls->rtol = opt->rtol;
  dls->atol = opt->atol;
/*--------------------------------------------------------------*/
  if((opt->itask < 2) || (opt->itask > 5))
    dls->itask = 1;
  else{
    dls->itask = opt->itask;
    if(dls->itask > 3)
      dls->rwork[0] = opt->t_critical;
  }
/*--------------------------------------------------------------*/
  dls->istate = 1;
/*--------------------------------------------------------------*/
  dls->jt = opt->jac_type;
  if((dls->jt == 4) || (dls->jt == 5)){
    dls->iwork[0] = opt->jac_lower_bandwidth;
    dls->iwork[1] = opt->jac_upper_bandwidth;
  }
/*--------------------------------------------------------------*/  
  return dls;
}

void dlsodar_workspace_reset(dlsodar_workspace *dls, const dlsodar_options *opt){
  int i, j, k;  
/*-----------------------------------------------------------------------*/
  int neq = opt->NEQ, ng = opt->NG;
/*-----------------------------------------------------------------------*/ 
  if(dls->ng < ng)
    dls->jroot = realloc(dls->jroot, sizeof(int) * ng);
  did_malloc_work(dls->jroot, "dlsodar_workspace_reset");

  for(i = 0; i < ng; i++)
    dls->jroot[i] = 0;
/*-----------------------------------------------------------------------*/
  dls->nerr = 0;
/*--------------------------------------------------------------*/
  int lrw_req;
  lrw_req = 22 + neq * max(16, neq + 9) + 3 * ng;

  if(dls->lrw < lrw_req)
    dls->rwork = realloc(dls->rwork, sizeof(double) * lrw_req);
  did_malloc_work(dls->rwork, "dlsodar_workspace_alloc");
  
  dls->lrw = lrw_req;
  
  for(i = 0; i < lrw_req; i++)
    dls->rwork[i] = 0.;
/*--------------------------------------------------------------*/
  double liw_req;
  liw_req = 20 + neq;

  if(dls->liw < liw_req)
    dls->iwork = realloc(dls->iwork, sizeof(int) * liw_req);  
  did_malloc_work(dls->iwork, "dlsodar_workspace_alloc");

  dls->liw = liw_req;

  for(i = 0; i < liw_req; i++)
    dls->iwork[i] = 0;
/*-----------------------------------------------------------------------*/
  dls->neq = neq;
  dls->ng = ng;
/*-----------------------------------------------------------------------*/  
  dls->iopt = 1;

  dls->rwork[4] = opt->init_step_size;
  dls->rwork[5] = opt->max_step_size;
  dls->rwork[6] = opt->min_step_size;

  dls->iwork[4] = opt->print_at_switch;
  dls->iwork[5] = opt->max_steps;
  dls->iwork[6] = opt->max_num_messages;
  dls->iwork[7] = opt->max_order_non_stiff;
  dls->iwork[8] = opt->max_order_stiff;
/*--------------------------------------------------------------*/
  if((opt->itol < 2) || (opt->itol > 4)){
    dls->itol = 1;//Error bounds are same for all dimensions.
  }
  else
    dls->itol = opt->itol;
  
  dls->rtol = opt->rtol;
  dls->atol = opt->atol;
/*--------------------------------------------------------------*/
  if((opt->itask < 2) || (opt->itask > 5))
    dls->itask = 1;
  else{
    dls->itask = opt->itask;
    if(dls->itask > 3)
      dls->rwork[0] = opt->t_critical;
  }
/*--------------------------------------------------------------*/
  dls->istate = 1;
/*--------------------------------------------------------------*/
  dls->jt = opt->jac_type;
  if((dls->jt == 4) || (dls->jt == 5)){
    dls->iwork[0] = opt->jac_lower_bandwidth;
    dls->iwork[1] = opt->jac_upper_bandwidth;
  }
/*--------------------------------------------------------------*/
  dls->nerr = 0;
  
  return;
}

void dlsodar_workspace_free(dlsodar_workspace *dls){
  free(dls->rwork);
  free(dls->iwork);
  free(dls->jroot);
  free(dls);
}

void dlsodar_integrate(double t,
		       double *t0, double *q,
		       odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func,
		       void *data, dlsodar_workspace *dls){

  /*Not ANSI, don't know if this is thread safe.*/
  void dlsodar_field_compat(const int *neq, const double *t_, const double *y, double *ydot){
    f_func(ydot, *t_, y, data);
    return;
  }
  
  void dlsodar_jacobian_compat(const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd){
    j_func(pd, *t_, y, data);
    return;
  }

  void dlsodar_constraint_compat(const int *neq, const double *t_, const double *y, const int *ng, double *gout){
    c_func(gout, *t_, y, data);
  }
  
  dlsodar_(dlsodar_field_compat,
	   &dls->neq, q, t0, &t,
	   &dls->itol, dls->rtol, dls->atol,
	   &dls->itask, &dls->istate, &dls->iopt,
	   dls->rwork, &dls->lrw, dls->iwork, &dls->liw,
	   dlsodar_jacobian_compat, &dls->jt,
	   dlsodar_constraint_compat, &dls->ng,
	   dls->jroot);
  
  if(dls->istate < 0){
    fprintf(stderr, "DLSODAR:\n Warning: ISTATE = %d\n", dls->istate);
    dls->nerr++;
  }
  return;
}
