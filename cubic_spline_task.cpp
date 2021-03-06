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
static double      start_time = 0.0;
static SL_DJstate  target[N_DOFS+1];
static double      delta_t = 0.01;
static double      duration = 1.0;
static double      time_to_go;


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


/*****************************************************************************
 ******************************************************************************
 Function Name    : add_cubic_spline_task
 Date        : Feb 1999
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
    
    addTask("Cubic Spline Task", init_cubic_spline_task,
            run_cubic_spline_task, change_cubic_spline_task);
    
}

/*****************************************************************************
 ******************************************************************************
 Function Name    : init_cubic_spline_task
 Date        : Dec. 1997
 
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
    static int firsttime = TRUE;
    
    if (firsttime){
        firsttime = FALSE;
    }
    
    // prepare going to the default posture
    bzero((char *)&(target[1]),N_DOFS*sizeof(target[1]));
    for (i=1; i<=N_DOFS; i++)
        target[i] = joint_default_state[i];
    
    // go to the target using inverse dynamics (ID)
    if (!go_target_wait_ID(target))
        return FALSE;
    
    // re-use the variable target for our min-jerk movement: only the right arm moves
    target[R_SFE].th += 0.4;
    target[R_SAA].th -= 0.4;
    target[R_EB].th  += 0.5;
    
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
 Function Name    : run_cubic_spline_task
 Date        : Dec. 1997
 
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
    if (time_to_go <= 0)
        freeze();
    
    return TRUE;
}

/*****************************************************************************
 ******************************************************************************
 Function Name    : change_cubic_spline_task
 Date        : Dec. 1997
 
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
