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

#ifndef __C_ODEPACK_DLSODAR__
#define __C_ODEPACK_DLSODAR__

/*---------------JT----------------------*/
#define USR_FULL_JAC 1
#define INT_FULL_JAC 2

#define USR_BAND_JAC 4
#define INT_BAND_JAC 5
/*-----------------------------------------------------------------------*/

typedef struct{
  int NEQ, NG;
  
  int jac_type;
  /* i - lower_bandwidth <= j <= i + upper_bandwidth */
  int jac_lower_bandwidth, jac_upper_bandwidth;

  int itol;
  double *rtol, *atol;
/*Quoting from opkdmain.f
C              ITOL    RTOL      ATOL      EWT(i)
C              ----    ------    ------    -----------------------------
C              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
C              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
C              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
C              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
*/   

  int itask;
  double t_critical;
/*Quoting from opkdmain.f
C     ITASK    An index specifying the task to be performed.  Input
C              only.  ITASK has the following values and meanings:
C              1   Normal computation of output values of y(t) at
C                  t = TOUT (by overshooting and interpolating).
C              2   Take one step only and return.
C              3   Stop at the first internal mesh point at or beyond
C                  t = TOUT and return.
C              4   Normal computation of output values of y(t) at
C                  t = TOUT but without overshooting t = TCRIT.  TCRIT
C                  must be input as RWORK(1).  TCRIT may be equal to or
C                  beyond TOUT, but not behind it in the direction of
C                  integration.  This option is useful if the problem
C                  has a singularity at or beyond t = TCRIT.
C              5   Take one step, without passing TCRIT, and return.
C                  TCRIT must be input as RWORK(1).
C
C              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C              (within roundoff), it will return T = TCRIT (exactly) to
C              indicate this (unless ITASK = 4 and TOUT comes before
C              TCRIT, in which case answers at T = TOUT are returned
C              first).
*/

  double init_step_size, max_step_size, min_step_size;
  int max_steps,\
    max_order_stiff, max_order_non_stiff,\
    print_at_switch, max_num_messages; 
/*
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C IXPR    IWORK(5)  flag to generate extra printing at method switches.
C                   IXPR = 0 means no extra printing (the default).
C                   IXPR = 1 means print data on each switch.
C                   T, H, and NST will be printed on the same logical
C                   unit as used for error messages.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
C                   (Adams) method.  The default value is 12.
C                   If MXORDN exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDN is held constant during the problem.
C
C MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
C                   (BDF) method.  The default value is 5.
C                   If MXORDS exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDS is held constant during the problem.
  */
} dlsodar_options;

typedef struct{
  int iopt, istate, itask, itol, *iwork, liw, lrw, neq, ng, nerr, *jroot, jt;
  double *atol, *rtol, *rwork;
} dlsodar_workspace;

dlsodar_workspace * dlsodar_workspace_alloc(const dlsodar_options *opt);

void dlsodar_workspace_reset(dlsodar_workspace *dls, const dlsodar_options *opt);

void dlsodar_workspace_free(dlsodar_workspace *dls);

void dlsodar_integrate(double t,
		       double *t0, double *q,
		       odepack_field_func f_func, odepack_jacobian_func j_func, odepack_root_func c_func,
		       void *data, dlsodar_workspace *dls);
#endif