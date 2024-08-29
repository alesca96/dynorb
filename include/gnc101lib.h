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
typedef void *odeSysFun(const double in_t, const double in_yy[], void *in_params, double out_dyydt[]);

/* Runge-Kutta 1 to 4 */
void gnc_rk1to4(odeSysFun in_sys, const double in_t0, const double in_t1, double in_h, const int in_rk_order, const double in_yy0[], double out_tt[], double out_yyt[])
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
        if (a == NULL || b == NULL || c == NULL)
        {
            fprintf(stderr, "Memory allocation failed.\n");
        }
        a[0] = 0.0;
        b[0] = 0.0;
        c[0] = 1.0;
        break;

    case 2: // RK2 (Heun)
        n_stages = 2;
        a = (double *)malloc(2 * sizeof(double));     // Note: a is [1 x 2]
        b = (double *)malloc(2 * 1 * sizeof(double)); // Note: b is [2 x 1]
        c = (double *)malloc(2 * sizeof(double));     // Note: c is [1 x 2]
        if (a == NULL || b == NULL || c == NULL)
        {
            fprintf(stderr, "Memory allocation failed.\n");
            break;
        }
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
        if (a == NULL || b == NULL || c == NULL)
        {
            fprintf(stderr, "Memory allocation failed.\n");
            break;
        }
        a[0] = 0.0;
        a[1] = 0.5;
        a[2] = 1.0;
        b[0] = 0.0;
        b[1] = 0.5;
        b[2] = -1.0;
        b[3] = 0.0;
        b[4] = 0.0;
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
        if (a == NULL || b == NULL || c == NULL)
        {
            fprintf(stderr, "Memory allocation failed.\n");
            break;
        }
        a[0] = 0.0;
        a[1] = 0.5;
        a[2] = 0.5;
        a[3] = 1.0;
        b[0] = 0.0;
        b[1] = 0.5;
        b[2] = 0.0;
        b[3] = 0.0;
        b[4] = 0.0;
        b[5] = 0.0;
        b[6] = 0.5;
        b[7] = 0.0;
        b[8] = 0.0;
        b[9] = 0.0;
        b[10] = 0.0;
        b[11] = 1.0;
        c[0] = 1.0 / 6;
        c[1] = 1.0 / 3;
        c[2] = 1.0 / 3;
        c[3] = 1.0 / 6;
        break;

    default: // Error
        printf("Error: The user-provided in_rk_order = %d is not valid. Please provide an order within the range [1, 4].\n", in_rk_order);
        break; // Exit if invalid input
    }

    // Allocate Memory for Current State and Copy in it Initial One:
    int sys_size = sizeof(in_yy0) / sizeof(in_yy0[0]);
    double *yy = (double *)malloc(sys_size * sizeof(double));
    double *yy_ti = (double *)malloc(sys_size * sizeof(double));
    memcpy(yy, in_yy0, sys_size * sizeof(double));

    // Allocate Memory for derivatives:
    double *dyydt = (double *)malloc(sys_size * sizeof(double));

    // Integration Time Instant:
    double t = in_t0;
    double ti = in_t0;

    // Allocate Memory for inner states:
    double *yy_inner = (double *)malloc(sys_size * sizeof(double));

    // Numerical Integration:
    while (t < in_t1)
    {
        // Current Time Inst. and State:
        ti = t;
        memcpy(yy_ti, yy, sys_size * sizeof(double));

        // Evaluate Time Derivatives at 'n_stages' points in [in_t0, in_t1]
        for (int i = 0; i < n_stages; ++i)
        {
            double t_inner = ti + a[i] * in_h;
            memcpy(yy_inner, yy_ti, sys_size * sizeof(double));
            for (int j = 0; j < (i - 1); ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    // TODO: yy_inner[k] += in_h * b[i + j] * f[:,j];
                }
                // TODO: in_sys(t_inner, yy_inner, NULL, f[:, i]);
            }
        }
        // Take min(h, tf-t)
        if (in_h > (in_t1 - t))
        {
            in_h = (in_t1 - t);
        }
        // Update t:
        t = t + in_h;
        // New State:
        // TODO: yy = yy_ti + in_h*f*c
        // out_tt = append(out_tt, t)
        // out_tt = append(out_yyt, yy)
    }

    // TODO:
    //      * matrxi_copy to append files
    //      *
    //
    // Free Memory:
    free(yy_inner);
    free(dyydt);
    free(yy_ti);
    free(yy);
    free(c);
    free(b);
    free(a);
}

#endif // GNC101LIB_C_

#ifdef GNCLIB_IMPLEMENTATION

#endif // GNCLIB_IMPLEMENTATION