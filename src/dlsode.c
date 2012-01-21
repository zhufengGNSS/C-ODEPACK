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

dlsode_workspace * dlsode_workspace_alloc(const dlsode_options *opt){
  dlsode_workspace *dls;
  dls = (dlsode_workspace *) calloc(1, sizeof(dlsode_workspace));
  did_malloc_work(dls, "dlsode_workspace_alloc");
  
  dls->neq = opt->NEQ;
  dls->mf = opt->step_method * 10 + opt->iter_method;

/*--------------------------------------------------------------*/  
  int max_ord, lwm;

  if(opt->step_method == ADAMS_IMPLICIT){
    if(opt->max_order == 0)
      max_ord = 12;
    else
      max_ord = min(opt->max_order, 12);
  }
  else if(opt->step_method == BDF){
    if(opt->max_order == 0)
      max_ord = 5;
    else
      max_ord = min(opt->max_order, 5);
  }
  else{
    fprintf(stderr, "dlsode_workspace_alloc: Unknown step method. METH = %d\nQuitting!\n", opt->step_method);
    exit(2);
  }

  if(opt->iter_method == 0){
    lwm = 0;
    dls->liw = 20;
  }
  else if((opt->iter_method == 1) || (opt->iter_method == 2)){
    lwm = opt->NEQ * opt->NEQ + 2;
    dls->liw = 20 + opt->NEQ;
  }
  else if(opt->iter_method == 3){
    lwm = opt->NEQ + 2;
    dls->liw = 20;
  }
  else if((opt->iter_method == 4) || (opt->iter_method == 5)){
    lwm = (2 * opt->jac_lower_bandwidth + opt->jac_upper_bandwidth + 1) * opt->NEQ + 2;
    dls->liw = 20 + opt->NEQ;
  }
  else{
    fprintf(stderr, "dlsode_workspace_alloc: Unknown iteration method. MITER = %d\nQuitting!\n", opt->iter_method);
    exit(2);
  }

  dls->lrw = 20 + opt->NEQ * (max_ord + 1) + 3 * opt->NEQ + lwm;
/*--------------------------------------------------------------*/
  dls->iopt = 1;
    
  dls->rwork = calloc(dls->lrw, sizeof(double));
  did_malloc_work(dls->rwork, "dlsode_workspace_alloc");
    
  dls->iwork = calloc(dls->liw, sizeof(int));
  did_malloc_work(dls->iwork, "dlsode_workspace_alloc");

  dls->rwork[4] = opt->init_step_size;
  dls->rwork[5] = opt->max_step_size;
  dls->rwork[6] = opt->min_step_size;
  
  dls->iwork[4] = opt->max_order;
  dls->iwork[5] = opt->max_steps;
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
  if((opt->iter_method == 4) || (opt->iter_method == 5)){
    dls->iwork[0] = opt->jac_lower_bandwidth;
    dls->iwork[1] = opt->jac_upper_bandwidth;
  }
/*--------------------------------------------------------------*/    
  return dls;
}

void dlsode_workspace_reset(dlsode_workspace *dls, const dlsode_options *opt){
  int i, j, k;
  
/*-----------------------------------------------------------------------*/  
  dls->neq = opt->NEQ;
  dls->mf = opt->step_method * 10 + opt->iter_method;
/*--------------------------------------------------------------*/  
  int max_ord, lwm, liw_req, lrw_req;

  if(opt->step_method == ADAMS_IMPLICIT){
    if(opt->max_order == 0)
      max_ord = 12;
    else
      max_ord = min(opt->max_order, 12);
  }
  else if(opt->step_method == BDF){
    if(opt->max_order == 0)
      max_ord = 5;
    else
      max_ord = min(opt->max_order, 5);
  }
  else{
    fprintf(stderr, "dlsode_workspace_reset: Unknown step method. METH = %d\nQuitting!\n", opt->step_method);
    exit(2);
  }

  if(opt->iter_method == 0){
    lwm = 0;
    liw_req = 20;
  }
  else if((opt->iter_method == 1) || (opt->iter_method == 2)){
    lwm = opt->NEQ * opt->NEQ + 2;
    liw_req = 20 + opt->NEQ;
  }
  else if(opt->iter_method == 3){
    lwm = opt->NEQ + 2;
    liw_req = 20;
  }
  else if((opt->iter_method == 4) || (opt->iter_method == 5)){
    lwm = (2 * opt->jac_lower_bandwidth + opt->jac_upper_bandwidth + 1) * opt->NEQ + 2;
    liw_req = 20 + opt->NEQ;
  }
  else{
    fprintf(stderr, "dlsode_workspace_reset: Unknown iteration method. MITER = %d\nQuitting!\n", opt->iter_method);
    exit(2);
  }

  lrw_req = 20 + opt->NEQ * (max_ord + 1) + 3 * opt->NEQ + lwm;
/*-----------------------------------------------------------------------*/
  if(lrw_req > dls->lrw){
    dls->rwork = realloc(dls->rwork, lrw_req * sizeof(double));
    did_malloc_work(dls->rwork, "dlsode_workspace_reset");
    dls->lrw = lrw_req;
  }

  for(i = 0; i < dls->lrw; i++)
    dls->rwork[i] = 0.;
/*-----------------------------------------------------------------------*/
  if(liw_req > dls->liw){
    dls->iwork = realloc(dls->iwork, liw_req * sizeof(int));
    did_malloc_work(dls->iwork, "dlsode_workspace_reset");
    dls->liw = liw_req;
  }

  for(i = 0; i < dls->liw; i++)
    dls->iwork[i] = 0.; 
/*-----------------------------------------------------------------------*/
  
  dls->iopt = 1;
  
  dls->rwork[4] = opt->init_step_size;
  dls->rwork[5] = opt->max_step_size;
  dls->rwork[6] = opt->min_step_size;
  
  dls->iwork[4] = opt->max_order;
  dls->iwork[5] = opt->max_steps;
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
  if((opt->iter_method == 4) || (opt->iter_method == 5)){
    dls->iwork[0] = opt->jac_lower_bandwidth;
    dls->iwork[1] = opt->jac_upper_bandwidth;
  }
/*--------------------------------------------------------------*/
  dls->nerr = 0;
  
  return;
}

void dlsode_workspace_free(dlsode_workspace *dls){
  free(dls->rwork);
  free(dls->iwork);
  free(dls);
}

void dlsode_integrate(double t,
		      double *t0, double *q,
		      odepack_field_func f_func, odepack_jacobian_func j_func,
		      void *data, dlsode_workspace *dls){

  /*Not ANSI, don't know if this is thread safe.*/
  void dlsode_field_compat(const int *neq, const double *t_, const double *y, double *ydot){
    f_func(ydot, *t_, y, data);
    return;
  }

  void dlsode_jacobian_compat(const int *neq, const double *t_, const double *y, const int *ml, const int *mu, double *pd, const int *nrowpd){
    j_func(pd, *t_, q, data);
    return;
  }

  dlsode_(dlsode_field_compat,
	  &dls->neq, q, t0, &t,
	  &dls->itol, dls->rtol, dls->atol,
	  &dls->itask, &dls->istate, &dls->iopt,
	  dls->rwork, &dls->lrw, dls->iwork, &dls->liw,
	  dlsode_jacobian_compat,
	  &dls->mf);
  
  if(dls->istate < 0){
    fprintf(stderr, "DLSODE:\n Warning: ISTATE = %d\n", dls->istate);
    dls->nerr++;
  }
  return;
}
