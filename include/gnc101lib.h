/* ****************************************************
 * **************** TODO: add description *************
 * **************************************************** */

/* ASSUMPTIONS:
 * 1. Code is written using COLUMN MAJOR convention (see CBLAS).

*/

#ifndef GNC101LIB_C_
#define GNC101LIB_C_

/* Includes: */
// TODO: #include <cblas.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For memcpy

/* Generic ODE Function: */
typedef void(_gnc_odeFun)(const double in_t, const double *in_yy, const void *in_params, double *out_dyydt);

/* ODE System Structure: */
typedef struct
{
    const _gnc_odeFun *odeFunction; // Pointer to the ODE function
    const void *params;             // Pointer to the parameters for the ODE function
    const double *yy0;              // Pointer to the initial state array
    const double t0;                // Initial time
    const double t1;                // Final time
    const int sys_size;             // Size of the system (number of equations)

} _gnc_odeSys;

/* **************************************************** */
void _gnc_rk1to4(_gnc_odeSys *in_sys, const int in_rk_order, const double in_h, const int in_n_steps, double *out_tt, double *out_yyt);
/*
 * Function for performing Runge-Kutta numerical integration
 * for ODE systems of various orders (RK1 to RK4).
 *
 * The input parameters:
 * - in_sys: Pointer to the structure defining the ODE system.
 * - in_rk_order: Order of the Runge-Kutta method (1, 2, 3, or 4).
 * - in_h: Time step size.
 * - out_tt: Output Pointer to array to store time points.
 * - out_yyt: Output Pointer to array to store the state vectors at each time point.
 *
 * This function uses a column-major convention.
 * **************************************************** */

#endif // GNC101LIB_C_

/*
 * **************************************************** *
 * ************* FUNCTIONS IMPLEMENTATION ************* *
 * **************************************************** *
 */

#ifdef GNCLIB_IMPLEMENTATION

void _gnc_rk1to4(_gnc_odeSys *in_sys, const int in_rk_order, const double in_h, const int in_n_steps, double *out_tt, double *out_yyt)
{
    // Open up _gnc_odeSys [TODO: remove this]
    _gnc_odeFun *odeFunction = in_sys->odeFunction; // Pointer to _gnc_odeFun
    const void *params = in_sys->params;            // Pointer to params
    const int sys_size = in_sys->sys_size;          // Size of System
    const double *yy0 = in_sys->yy0;                // Pointer to Initial Conditions
    double t0 = in_sys->t0;                         // Initial time
    double t1 = in_sys->t1;                         // Final time

    // Declare n_stages and pointers for a, b, c
    int n_stages = 0;
    double *a = NULL;
    double *b = NULL;
    double *c = NULL;

    // Selection of RK1, RK2, ..., RK4:
    switch (in_rk_order)
    {
    case 1: // RK1 (Euler)
        n_stages = 1;
        a = (double *)malloc(1 * sizeof(double));     // Note: a is [1 x 1]
        b = (double *)malloc(1 * 1 * sizeof(double)); // Note: b is [1 x 1]
        c = (double *)malloc(1 * sizeof(double));     // Note: c is [1 x 1]

        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = 1.0;
        break;

    case 2: // RK2 (Heun)
        n_stages = 2;
        a = (double *)malloc(2 * sizeof(double));     // Note: a is [1 x 2]
        b = (double *)malloc(2 * 1 * sizeof(double)); // Note: b is [2 x 1]
        c = (double *)malloc(2 * sizeof(double));     // Note: c is [2 x 1]
        a[0] = 0.0;
        a[1] = 1.0;
        b[0] = 0.0;
        b[1] = 1.0;
        c[0] = 0.5;
        c[1] = 0.5;
        break;

    case 3: // RK3
        n_stages = 3;
        a = (double *)malloc(3 * sizeof(double));     // Note: a is [1 x 3]
        b = (double *)malloc(3 * 2 * sizeof(double)); // Note: b is [3 x 2]
        c = (double *)malloc(3 * sizeof(double));     // Note: c is [3 x 1]
        a[0] = 0.0;
        a[1] = 0.5;
        a[2] = 1.0;
        b[0] = 0.0;
        b[1] = 0.0;
        b[2] = 0.5;
        b[3] = 0.0;
        b[4] = -1.0;
        b[5] = 2.0;
        c[0] = 1.0 / 6;
        c[1] = 2.0 / 3;
        c[2] = 1.0 / 6;
        break;

    case 4: // RK4
        n_stages = 4;
        a = (double *)malloc(4 * sizeof(double));     // Note: a is [1 x 4]
        b = (double *)malloc(4 * 3 * sizeof(double)); // Note: b is [4 x 3]
        c = (double *)malloc(4 * sizeof(double));     // Note: c is [4 x 1]

        a[0] = 0.0;
        a[1] = 0.5;
        a[2] = 0.5;
        a[3] = 1.0;
        b[0] = 0.0;
        b[1] = 0.0;
        b[2] = 0.0;
        b[3] = 0.5;
        b[4] = 0.0;
        b[5] = 0.0;
        b[6] = 0.0;
        b[7] = 0.5;
        b[8] = 0.0;
        b[9] = 0.0;
        b[10] = 0.0;
        b[11] = 1.0;
        c[0] = 1.0 / 6;
        c[1] = 1.0 / 3;
        c[2] = 1.0 / 3;
        c[3] = 1.0 / 6;
        break;

    default:
        printf("Error: The user-provided in_rk_order = %d is not valid. Please provide an order within the range [1, 4].\n", in_rk_order);
        return;
    }

    // Allocate Memory for Current State, Inner State:
    double *yy = (double *)malloc(sys_size * sizeof(double));
    double *yy_inner = (double *)malloc(sys_size * sizeof(double));
    memcpy(yy, yy0, sys_size * sizeof(double)); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    double *ff = (double *)malloc(sys_size * n_stages * sizeof(double)); // (ff = dyy/dt)

    // Integration Time Instant:
    double t = t0;

    // Numerical Integration:
    for (int step = 0; step < in_n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            double t_inner = t + a[i] * in_h;
            memcpy(yy_inner, yy, sys_size * sizeof(double));
            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += in_h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, params, &ff[i * sys_size]);
        }

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += in_h * c[i] * ff[i * sys_size + k];
            }
        }

        // Update Integration Time:
        t += in_h;

        // Store results:
        if (t > t1)
        {
            memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(double));
            out_tt[step] = t;
            printf("_gnc_rk1to4 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(double));
        out_tt[step] = t;
    }

    // Free Memory:
    // printf("_gnc_rk1to4: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    free(c);
    free(b);
    free(a);
    // printf("_gnc_rk1to4: Done Freeing Memory:\n");
}

#endif // GNCLIB_IMPLEMENTATION