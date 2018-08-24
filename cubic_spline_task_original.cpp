/*============================================================================
==============================================================================
                      
                              cubic_spline_task.c
 
==============================================================================
Remarks:

      sekeleton to create the sample task

============================================================================*/

// system headers
#include "SL_system_headers.h"

/* SL includes */
#include "SL.h"
#include "SL_user.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_kinematics.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_man.h"

// defines

// local variables
static double      time_step;
//static double      start_time = 0.0;
static SL_DJstate  target[N_DOFS+1];
//static double      delta_t = 0.01;
//static double      duration = 1.0;
//static double      time_to_go;
//static int         flag=0;
static double     *cart;
static SL_Cstate  *ctarget;
static SL_Cstate  *cnext;
static int        *cstatus;
static SL_DJstate *target;
static int         firsttime = TRUE;
static double      movement_time = 1.0;
static double      tau;

// global functions 
extern "C" void
add_cubic_spline_task( void );

// local functions
static int  init_cubic_spline_task(void);
static int  run_cubic_spline_task(void);
static int  change_cubic_spline_task(void);

static int 
cubic_spline_next_step (double x,double xd, double xdd, double t, double td, double tdd,
			double t_togo, double dt,
			double *x_next, double *xd_next, double *xdd_next);

static void init_vars(void);
static int  calculate_min_jerk_next_step (SL_Cstate *curr_state,
                                          SL_Cstate *des_state,
                                          double tau,
                                          double delta_t,
                                          SL_Cstate *next_states);


/*****************************************************************************
******************************************************************************
Function Name	: add_cubic_spline_task
Date		: Feb 1999
Remarks:
 
adds the task to the task menu

******************************************************************************
Paramters:  (i/o = input/output)

none

*****************************************************************************/
void
add_cubic_spline_task( void )
{
  int i, j;
  static int firsttime = TRUE;
    if (firsttime) {
        firsttime = FALSE;
        
        cart    = my_vector(1,n_endeffs*6);
        ctarget = (SL_Cstate *) my_calloc(n_endeffs+1,sizeof(SL_Cstate),MY_STOP);
        cnext   = (SL_Cstate *) my_calloc(n_endeffs+1,sizeof(SL_Cstate),MY_STOP);
        cstatus = my_ivector(1,n_endeffs*6);
        target  = (SL_DJstate *)my_calloc(n_dofs+1,sizeof(SL_DJstate),MY_STOP);

  
  addTask("Cubic Spline Task", init_cubic_spline_task, 
	  run_cubic_spline_task, change_cubic_spline_task);
    }

}    

/*****************************************************************************
******************************************************************************
  Function Name	: init_cubic_spline_task
  Date		: Dec. 1997

  Remarks:

  initialization for task

******************************************************************************
  Paramters:  (i/o = input/output)

       none

 *****************************************************************************/
static int 
init_cubic_spline_task(void)
{
  int j, i;
  int ans;
    char   string[100];
    double max_range=0;
    double aux;
    int    flag = FALSE;
    int    iaux;
    /*
  static int firsttime = TRUE;
  
  if (firsttime){
    firsttime = FALSE;
  }
*/
    /* check whether any other task is running */
    if (strcmp(current_task_name,NO_TASK) != 0) {
        printf("Goto task can only be run if no other task is running!\n");
        return FALSE;
    }
    
    /* initialize some variables */
    init_vars();
    time_step = 1./(double)task_servo_rate;
  // prepare going to the default posture
  bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
  for (i=1; i<=N_DOFS; i++)
    target[i] = joint_default_state[i];
  // go to the target using inverse dynamics (ID)
  if (!go_target_wait_ID(target)) 
    return FALSE;
// starting point for left hand
    target[L_SFE].th = 0;
    //target[L_SAA].th = -1 ;
    target[L_EB].th  = 1.44 ;
    target[L_HR].th  = 0;
//  // re-use the variable target for our min-jerk movement: only the right arm moves
//  target[L_SFE].th += 0.4;
//  target[L_SAA].th -= 0.4;
//  target[L_EB].th  += 0.5;


    
  // ready to go
  ans = 999;
  while (ans == 999) {
    if (!get_int("Enter 1 to start or anthing else to abort ...",ans,&ans))
      return FALSE;
  }
  
  // only go when user really types the right thing
  if (ans != 1) 
    return FALSE;

  start_time = task_servo_time;
  printf("start time = %.3f, task_servo_time = %.3f\n", 
	 start_time, task_servo_time);

  // start data collection
  scd();

  // time to go
  time_to_go = duration;

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: run_cubic_spline_task
  Date		: Dec. 1997

  Remarks:

  run the task from the task servo: REAL TIME requirements!

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
run_cubic_spline_task(void)
{
  int j, i;

  double task_time;

  // ******************************************
  // NOTE: all array indices start with 1 in SL
  // ******************************************

  task_time = task_servo_time - start_time;

  // compute the update for the desired states
  for (i=1; i<=N_DOFS; ++i) {
    cubic_spline_next_step(joint_des_state[i].th,
			   joint_des_state[i].thd,
			   joint_des_state[i].thdd,
			   target[i].th,
			   target[i].thd,
			   target[i].thdd,
			   time_to_go,
			   delta_t,
			   &(joint_des_state[i].th),
			   &(joint_des_state[i].thd),
			   &(joint_des_state[i].thdd));
  }

  // compute inverse dynamics torques
  SL_InvDynNE(joint_state,joint_des_state,endeff,&base_state,&base_orient);
  
  // decrement time to go
  time_to_go -= delta_t;
    if (time_to_go <= 0.02){
       //freeze();
        time_to_go = duration;
        if (flag < 4) {
        switch (flag){
            case 0:
                target[L_SFE].th = 0;
                target[L_SAA].th = -1.03;
                target[L_EB].th  = 1.43;
                target[L_HR].th  = 0;
                break;
            case 1:
                target[L_SFE].th = 0;
                target[L_SAA].th = -0.89;
                target[L_EB].th  = 1.11 ;
                target[L_HR].th  = 0;
                break;
            case 2:
                target[L_SFE].th = 0;
                target[L_SAA].th = -0.68 ;
                target[L_EB].th  = 0.97;
                target[L_HR].th  = 0;
                break;
            case 3:
                target[L_SFE].th = 0;
                target[L_SAA].th = -0.62;
                target[L_EB].th  = 1.21 ;
                target[L_HR].th  = 0;
                break;
        }
            flag++;
        } else {
            flag = 0;
        }
    }
    
  

  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: change_cubic_spline_task
  Date		: Dec. 1997

  Remarks:

  changes the task parameters

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
change_cubic_spline_task(void)
{
  int    ivar;
  double dvar;

  get_int("This is how to enter an integer variable",ivar,&ivar);
  get_double("This is how to enter a double variable",dvar,&dvar);

  return TRUE;

}


/*!*****************************************************************************
 *******************************************************************************
\note  cubic_spline_next_step
\date  April 2014
   
\remarks 

Given the time to go, the current state is updated to the next state
using min jerk splines

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]          x,xd,xdd : the current state, vel, acceleration
 \param[in]          t,td,tdd : the target state, vel, acceleration
 \param[in]          t_togo   : time to go until target is reached
 \param[in]          dt       : time increment
 \param[in]          x_next,xd_next,xdd_next : the next state after dt

 ******************************************************************************/
static int 
cubic_spline_next_step (double x,double xd, double xdd, double t, double td, double tdd,
			double t_togo, double dt,
			double *x_next, double *xd_next, double *xdd_next)

{
    double c0, c1, c2, c3;
  // your code goes here ...

        c0 = x;
        c1 = xd;
        c2 = (3 * t - 3 * x - 2 * xd * t_togo - td * t_togo) / (pow(t_togo ,2));
        c3 = (td * t_togo - 2 * t + 2 * x + xd * t_togo) / (pow(t_togo, 3));
        *x_next  = c0 + c1 * dt + c2 * pow(dt, 2) + c3 * pow(dt, 3);
        *xd_next = c1 + c2 * 2 * dt + c3 * 3 * pow(dt, 2);
        *xdd_next = c2 * 2 + c3 * dt * 6;
        x = *x_next;
        xd = *xd_next;
        xdd = *xdd_next;
  return TRUE;
}
/*!*****************************************************************************
 *******************************************************************************
 \note  calculate_min_jerk_next_step
 \date  August 1994
 
 \remarks
 
 given a current cart state and a target cart
 state as well as movement duration, the next increment
 for the cart states is calculated. Note that this
 is done in only in cartesian dimensions with active status.
 NOTE that this function requires velocity and accelerations
 as input as well!!!
 
 *******************************************************************************
 Function Parameters: [in]=input,[out]=output
 
 \param[in]     curr_states: the current state
 \param[in]     des_states : the desired state
 \param[in]     tau        : the desired movement duration until the goal is
 reached.
 \param[in]     dt         : at which delta time is the next_states from now
 \param[out]    next_states: the next state after dt
 
 ******************************************************************************/
static int
calculate_min_jerk_next_step (SL_Cstate *curr_state,
                              SL_Cstate *des_state,
                              double tau,
                              double delta_t,
                              SL_Cstate *next_state)

{
    double t1,t2,t3,t4,t5;
    double tau1,tau2,tau3,tau4,tau5;
    int    i,j;
    
    if (delta_t > tau || tau < 1./(double)task_servo_rate || delta_t <= 0) {
        return FALSE;
    }
    
    t1 = delta_t;
    t2 = t1 * delta_t;
    t3 = t2 * delta_t;
    t4 = t3 * delta_t;
    t5 = t4 * delta_t;
    
    tau1 = tau;
    tau2 = tau1 * tau;
    tau3 = tau2 * tau;
    tau4 = tau3 * tau;
    tau5 = tau4 * tau;
    
    for (j=1; j<=n_endeffs; ++j) {
        for (i=1; i<=N_CART; ++i) {
            
            if (cstatus[(j-1)*6+i]) {
                
                /* calculate the constants */
                
                const double dist   = des_state[j].x[i] - curr_state[j].x[i];
                const double p1     = des_state[j].x[i];
                const double p0     = curr_state[j].x[i];
                const double a1t2   = des_state[j].xdd[i];
                const double a0t2   = curr_state[j].xdd[i];
                const double v1t1   = des_state[j].xd[i];
                const double v0t1   = curr_state[j].xd[i];
                
                const double c1 = 6.*dist/tau5 + (a1t2 - a0t2)/(2.*tau3) -
                3.*(v0t1 + v1t1)/tau4;
                const double c2 = -15.*dist/tau4 + (3.*a0t2 - 2.*a1t2)/(2.*tau2) +
                (8.*v0t1 + 7.*v1t1)/tau3;
                const double c3 = 10.*dist/tau3+ (a1t2 - 3.*a0t2)/(2.*tau) -
                (6.*v0t1 + 4.*v1t1)/tau2;
                const double c4 = curr_state[j].xdd[i]/2.;
                const double c5 = curr_state[j].xd[i];
                const double c6 = curr_state[j].x[i];
                
                next_state[j].x[i]   = c1*t5 + c2*t4 + c3*t3 + c4*t2 + c5*t1 + c6;
                next_state[j].xd[i]  = 5.*c1*t4 + 4*c2*t3 + 3*c3*t2 + 2*c4*t1 + c5;
                next_state[j].xdd[i] = 20.*c1*t3 + 12.*c2*t2 + 6.*c3*t1 + 2.*c4;
                
            }
        }
    }
    
    return TRUE;
}




