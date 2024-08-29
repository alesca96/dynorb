/* ****************************************************
 * **************** TODO: add description *************
 * **************************************************** */

/* ASSUMPTIONS:
 * 1. Code is written using COLUMN MAJOR convention (see CBLAS).

*/

#ifndef GNC101LIB_C_
#define GNC101LIB_C_

/* Includes: */
#include <stdlib.h>
#include <string.h> // For memcpy
#include <stdio.h>
#include <cblas.h>

/* Structures: */

/* Generic ODE system: */
typedef void (*odeSysFun)(const double in_t, const double in_yy[], void *in_params, double out_dyydt[]);

void gnc_rk1to4(odeSysFun in_sys, const double in_t0, const double in_t1, double in_h, const int in_rk_order, const double in_yy0[], int sys_size, double out_tt[], double out_yyt[]);

/* ****************************************************
* void gnc_rk1to4(odeSysFun in_sys, const double in_t0, const double in_t1,
                 double in_h, const int in_rk_order, const double in_yy0[],
                 double out_tt[], double out_yyt[])

 * Function for performing Runge-Kutta numerical integration
 * for ODE systems of various orders (RK1 to RK4).
 *
 * The input parameters:
 * - in_sys: Pointer to the function defining the ODE system.
 * - in_t0: Initial time.
 * - in_t1: Final time.
 * - in_h: Time step size.
 * - in_rk_order: Order of the Runge-Kutta method (1, 2, 3, or 4).
 * - in_yy0: Initial state vector at time in_t0.
 * - out_tt: Output array to store time points.
 * - out_yyt: Output array to store the state vectors at each time point.
 *
 * This function uses a column-major convention.
 * **************************************************** */

#endif // GNC101LIB_C_

#ifdef GNCLIB_IMPLEMENTATION

void gnc_rk1to4(odeSysFun in_sys, const double in_t0, const double in_t1, double in_h, const int in_rk_order, const double in_yy0[], int sys_size, double out_tt[], double out_yyt[])
{
    int n_stages = 0; // Declare n_stages outside of switch
    double *a = NULL; // Declare pointers for a, b, c
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

    // Allocate Memory for Current State and Copy Initial One:
    double *yy = (double *)malloc(sys_size * sizeof(double));
    double *yy_inner = (double *)malloc(sys_size * sizeof(double));
    memcpy(yy, in_yy0, sys_size * sizeof(double));

    // Allocate Memory for derivatives:
    double *ff = (double *)malloc(sys_size * n_stages * sizeof(double));

    // Integration Time Instant:
    double t = in_t0;
    int step = 0;

    // Store initial state:
    memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(double));
    out_tt[step] = t;
    step++;

    // Numerical Integration:
    while (t < in_t1)
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
            in_sys(t_inner, yy_inner, NULL, &ff[i * sys_size]);
        }

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += in_h * c[i] * ff[i * sys_size + k];
            }
        }

        // Update time and store results:
        t += in_h;
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(double));
        out_tt[step] = t;
        step++;
    }

    // Free Memory:
    free(yy_inner);
    free(ff);
    free(yy);
    free(c);
    free(b);
    free(a);
}

#endif // GNCLIB_IMPLEMENTATION