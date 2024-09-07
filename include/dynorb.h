/* ****************************************************
 * ****************** DYNORB HEADER *******************
 * **************************************************** */

/* ASSUMPTIONS:
 * 1. Code uses COLUMN MAJOR convention (similar to CBLAS).
 * 2. Strides other than 1 are not supported in this implementation.
 */

#ifndef DYNORB_H_
#define DYNORB_H_

/* Double-Float switch: */
#ifdef USE_DOUBLE
typedef double real;
#elif defined(USE_FLOAT)
typedef float real;
#else
#error "Either USE_DOUBLE or USE_FLOAT must be defined"
#endif

/* CBLAS switch: */
#ifdef USE_CBLAS
#include <cblas.h>
#ifdef USE_DOUBLE
#define _dynorb_raxpy cblas_daxpy /* Double precision */
#elif defined(USE_FLOAT)
#define _dynorb_raxpy cblas_saxpy /* Single precision */
#endif

#else

/* Fallback CBLAS-like enums: */
typedef enum CBLAS_LAYOUT
{
    CblasRowMajor = 101,
    CblasColMajor = 102
} CBLAS_LAYOUT;

typedef enum CBLAS_TRANSPOSE
{
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113
} CBLAS_TRANSPOSE;

typedef enum CBLAS_UPLO
{
    CblasUpper = 121,
    CblasLower = 122
} CBLAS_UPLO;

typedef enum CBLAS_DIAG
{
    CblasNonUnit = 131,
    CblasUnit = 132
} CBLAS_DIAG;

typedef enum CBLAS_SIDE
{
    CblasLeft = 141,
    CblasRight = 142
} CBLAS_SIDE;
#endif

/* Standard C library includes: */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* Fallback CBLAS-like functions: */
#ifndef USE_CBLAS
void stride_warning(void);
void _dynorb_raxpy(const int N, const real alpha, const real *X,
                   const int incX, real *Y, const int incY);
#endif

/* MACROS: */
#define EL(M, nrows, ncols, i, j) M[(j * (nrows)) + i]

/* Utility functions: */
real min(real a, real b);
real max(real a, real b);
void matcpy(real *dst_M, const real *src_M, const int start_idx_row, const int start_idx_col, const int end_idx_row, const int end_idx_col, const int src_ncols, const int dst_ncols);
void matprint(const real *in_M, const int nrows, const int ncols);

/* ODE Function Type: */
typedef void(_dynorb_odeFun)(const real in_t, const real *in_yy, const void *in_params, real *out_dyydt);

/* ODE System Structure Type: */
typedef struct
{
    _dynorb_odeFun *odeFunction; // Pointer to the ODE function
    const void *params;          // Pointer to the parameters for the ODE function
    const real *yy0;             // Pointer to the initial state array
    const real t0;               // Initial time
    const real t1;               // Final time
    const int sys_size;          // Size of the system (number of equations)
    real *tt;                    // Time Steps of the Solution
    real *yyt;                   // Solution Array

} _dynorb_odeSys;

/* Function Declaration: */

/**
 *
 * @brief Performs Euler (RK1) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Euler (Runge-Kutta order 1) method.
 *
 * @param[in] in_sys Pointer to the structure defining the ODE system.
 * @param[in] in_h Time step size for the integration.
 * @param[in] in_n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk1(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps);

/**
 *
 * @brief Performs Heun (RK2) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Heun (Runge-Kutta order 2) method.
 *
 * @param[in] in_sys Pointer to the structure defining the ODE system.
 * @param[in] in_h Time step size for the integration.
 * @param[in] in_n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk2(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps);

/**
 *
 * @brief Performs Heun (RK3) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Runge-Kutta order 3 method.
 *
 * @param[in] in_sys Pointer to the structure defining the ODE system.
 * @param[in] in_h Time step size for the integration.
 * @param[in] in_n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk3(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps);

/**
 *
 * @brief Performs Heun (RK4) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Runge-Kutta order 4 method.
 *
 * @param[in] in_sys Pointer to the structure defining the ODE system.
 * @param[in] in_h Time step size for the integration.
 * @param[in] in_n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk4(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps);

/**
 *
 * @brief Performs Heun Predictor-Corrector numerical integration for ODE systems.
 *
 * This function integrates an ODE system using Heun Predictor-Corrector method.
 *
 * @param[in] in_sys Pointer to the structure defining the ODE system.
 * @param[in] in_h Time step size for the integration.
 * @param[in] in_n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_heun_(_dynorb_odeSys *in_sys, real in_h, const int in_n_steps);

#endif // DYNORB_H_

/*
 * **************************************************** *
 * ************* DYNORN IMPLEMENTATION ************* *
 * **************************************************** *
 */

#ifdef DYNORB_IMPLEMENTATION

/* Fallback CBLAS-like functions: */
#ifndef USE_CBLAS /* fallback to my implementation-> */
void stride_warning(void)
{
    printf("\nWARNING: Stride (inc_) different from 1 not supported.\n");
    printf("         Stride set to default: 1.\n");
    printf("         To use CBLAS: '#define USE_CBLAS'.\n\n");
}

/* Vector Sum: */
void _dynorb_raxpy(const int N, const real alpha, const real *X,
                   const int incX, real *Y, const int incY)
{
    if (incY != 1 || incX != 1)
    {
        stride_warning();
    }

    for (int i = 0; i < N; ++i)
    {
        Y[i] += (alpha * X[i]);
    }
}
#endif

/* Utility functions: */
real min(real a, real b)
{
    return (a < b) ? a : b;
}

real max(real a, real b)
{
    return (a > b) ? a : b;
}

void matcpy(real *dst_M, const real *src_M, const int start_idx_row, const int start_idx_col, const int end_idx_row, const int end_idx_col, const int src_ncols, const int dst_ncols)
{
    // TODO: Make more Powrful by allowing copying in a portion of the dest.
    for (int i = start_idx_row; i <= end_idx_row; ++i)
    {
        for (int j = start_idx_col; j <= end_idx_col; ++j)
        {
            // Ensure the destination and source have the same number of columns for the operation
            if (i < 0 || j < 0 || i >= (end_idx_row + 1) || j >= (end_idx_col + 1))
            {
                fprintf(stderr, "Index out of bounds\n");
                return;
            }
            EL(dst_M, dst_ncols, end_idx_row - start_idx_row + 1, i - start_idx_row, j - start_idx_col) = EL(src_M, src_ncols, end_idx_row - start_idx_row + 1, i, j);
        }
    }
}

void matprint(const real *in_M, const int nrows, const int ncols)
{
    printf("[\n");
    for (int i = 0; i < nrows; ++i)
    {
        printf("    ");
        for (int j = 0; j < ncols; ++j)
        {
            printf("%.4f", EL(in_M, nrows, ncols, i, j));
            if (j < ncols - 1)
            {
                printf(", ");
            }
        }
        printf(";\n");
    }
    printf("]\n");
}
/* Core Library Functions: */
void _dynorb_rk1(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = in_sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = in_sys->params;               // Pointer to params
    const int sys_size = in_sys->sys_size;             // Size of System
    const real *yy0 = in_sys->yy0;                     // Pointer to Initial Conditions
    real t0 = in_sys->t0;                              // Initial time
    real t1 = in_sys->t1;                              // Final time
    real *out_tt = in_sys->tt;                         // Time steps of solution
    real *out_yyt = in_sys->yyt;                       // Solution Array

    //  RK1 (Euler)
    int n_stages = 1;
    real a = 0.0; // Note: a is [1 x 1]
    // real b = 1.0; // Note: b is [1 x 1]
    real c = 1.0; // Note: c is [1 x 1]

    // Allocate Memory for Current State, Inner State:
    real *yy = (real *)malloc(sys_size * sizeof(real));
    memcpy(yy, yy0, sys_size * sizeof(real)); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real *ff = (real *)malloc(sys_size * n_stages * sizeof(real)); // (ff = dyy/dt)

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < in_n_steps; ++step)
    {
        // Evaluate Time Derivatives in [t, t+h]
        real t_inner = t + a * in_h;
        odeFunction(t_inner, yy, params, ff);

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += in_h * c * ff[k];
            }
        }

        // Update Integration Time:
        t += in_h;

        // Store results:
        if (t > t1)
        {
            memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
            out_tt[step] = t;
            printf("_dynorb_rk1 : Breaking From Loop<t = %f [s]>", t);
        }
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        out_tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk1: Begin Freeing Memory:\n");
    free(ff);
    free(yy);
    // printf("_dynorb_rk1: Done Freeing Memory:\n");
}

void _dynorb_rk2(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = in_sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = in_sys->params;               // Pointer to params
    const int sys_size = in_sys->sys_size;             // Size of System
    const real *yy0 = in_sys->yy0;                     // Pointer to Initial Conditions
    real t0 = in_sys->t0;                              // Initial time
    real t1 = in_sys->t1;                              // Final time
    real *out_tt = in_sys->tt;                         // Time steps of solution
    real *out_yyt = in_sys->yyt;                       // Solution Array

    //  RK2 (Heun)
    int n_stages = 2;
    real a[2] = {0.0, 0.0}; // Note: a is [1 x 2]
    real b[2] = {0.0, 0.0}; // Note: b is [2 x 1]
    real c[2] = {0.0, 0.0}; // Note: c is [2 x 1]
    a[1] = 1.0;
    b[0] = 0.0;
    b[1] = 1.0;
    c[0] = 0.5;
    c[1] = 0.5;

    // Allocate Memory for Current State, Inner State:
    real *yy = (real *)malloc(sys_size * sizeof(real));
    real *yy_inner = (real *)malloc(sys_size * sizeof(real));
    memcpy(yy, yy0, sys_size * sizeof(real)); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real *ff = (real *)malloc(sys_size * n_stages * sizeof(real)); // (ff = dyy/dt)

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < in_n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * in_h;
            memcpy(yy_inner, yy, sys_size * sizeof(real));
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
            memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
            out_tt[step] = t;
            printf("_dynorb_rk2 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        out_tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk2: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    // printf("_dynorb_rk2: Done Freeing Memory:\n");
}

void _dynorb_rk3(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = in_sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = in_sys->params;               // Pointer to params
    const int sys_size = in_sys->sys_size;             // Size of System
    const real *yy0 = in_sys->yy0;                     // Pointer to Initial Conditions
    real t0 = in_sys->t0;                              // Initial time
    real t1 = in_sys->t1;                              // Final time
    real *out_tt = in_sys->tt;                         // Time steps of solution
    real *out_yyt = in_sys->yyt;                       // Solution Array

    //  RK3
    int n_stages = 3;
    real a[3];     // Note: a is [1 x 3]
    real b[3 * 2]; // Note: b is [3 x 2]
    real c[3];     // Note: c is [3 x 1]
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

    // Allocate Memory for Current State, Inner State:
    real *yy = (real *)malloc(sys_size * sizeof(real));
    real *yy_inner = (real *)malloc(sys_size * sizeof(real));
    memcpy(yy, yy0, sys_size * sizeof(real)); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real *ff = (real *)malloc(sys_size * n_stages * sizeof(real)); // (ff = dyy/dt)

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < in_n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * in_h;
            memcpy(yy_inner, yy, sys_size * sizeof(real));
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
            memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
            out_tt[step] = t;
            printf("_dynorb_rk3 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        out_tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk3: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    // printf("_dynorb_rk3: Done Freeing Memory:\n");
}

void _dynorb_rk4(_dynorb_odeSys *in_sys, const real in_h, const int in_n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = in_sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = in_sys->params;               // Pointer to params
    const int sys_size = in_sys->sys_size;             // Size of System
    const real *yy0 = in_sys->yy0;                     // Pointer to Initial Conditions
    real t0 = in_sys->t0;                              // Initial time
    real t1 = in_sys->t1;                              // Final time
    real *out_tt = in_sys->tt;                         // Time steps of solution
    real *out_yyt = in_sys->yyt;                       // Solution Array

    //  RK4 (Runge-Kutta)
    int n_stages = 4;
    real a[4];     // Note: a is [1 x 4]
    real b[4 * 3]; // Note: b is [4 x 3]
    real c[4];     // Note: c is [4 x 1]
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

    // Allocate Memory for Current State, Inner State:
    real *yy = (real *)malloc(sys_size * sizeof(real));
    real *yy_inner = (real *)malloc(sys_size * sizeof(real));
    memcpy(yy, yy0, sys_size * sizeof(real)); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real *ff = (real *)malloc(sys_size * n_stages * sizeof(real)); // (ff = dyy/dt)

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < in_n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * in_h;
            memcpy(yy_inner, yy, sys_size * sizeof(real));
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
            memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
            out_tt[step] = t;
            printf("_dynorb_rk4 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        out_tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk4: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    // printf("_dynorb_rk4: Done Freeing Memory:\n");
}

void _dynorb_heun_(_dynorb_odeSys *in_sys, real in_h, const int in_n_steps)
{
    // Open up _dynorb_odeSys
    _dynorb_odeFun *odeFunction = in_sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = in_sys->params;               // Pointer to params
    const int sys_size = in_sys->sys_size;             // Size of System
    const real *yy0 = in_sys->yy0;                     // Pointer to Initial Conditions
    real t0 = in_sys->t0;                              // Initial time
    real t1 = in_sys->t1;                              // Final time
    real *out_tt = in_sys->tt;                         // Time steps of solution
    real *out_yyt = in_sys->yyt;                       // Solution Array

    // Tolerance and Max number of steps:
    const real tol = 1.e-6;
    const int max_iter = 100;

    // Initialize Algorithm:
    real t = t0;
    real *yy = (real *)malloc(sys_size * sizeof(real));
    if (yy == NULL)
    {
        perror("Memory allocation failed for yy");
        exit(EXIT_FAILURE);
    }
    memcpy(yy, yy0, sys_size * sizeof(real));

    // Copy First State and Instant in Output:
    out_tt[0] = t0;
    memcpy(out_yyt, yy0, sys_size * sizeof(real));

    // Allocate Memory for state and derivatives at interval boundaries:
    real *yy1_ = (real *)malloc(sys_size * sizeof(real));
    real *yy2_ = (real *)malloc(sys_size * sizeof(real));
    real *ff1_ = (real *)malloc(sys_size * sizeof(real));
    real *ff2_ = (real *)malloc(sys_size * sizeof(real));
    real *yy2pred_ = (real *)malloc(sys_size * sizeof(real));
    real *ffavg_ = (real *)malloc(sys_size * sizeof(real));

    // Main Loop
    int step = 1;
    while (t < t1 && step < in_n_steps)
    {
        // Step Size:
        in_h = min(in_h, t1 - t);

        // Left Boundary of Interval:
        real t1_ = t;
        memcpy(yy1_, yy, sys_size * sizeof(real));
        odeFunction(t1_, yy1_, params, ff1_);

        // Compute yy2_
        for (int i = 0; i < sys_size; ++i)
        {
            yy2_[i] = yy1_[i] + ff1_[i] * in_h;
        }

        // Right Boundary of the Interval:
        real t2_ = t1_ + in_h;
        real err = tol + 1;
        int iter = 0;

        // Predictor-Corrector Loop
        while (err > tol && iter <= max_iter)
        {
            memcpy(yy2pred_, yy2_, sys_size * sizeof(real));
            odeFunction(t2_, yy2pred_, params, ff2_);

            // Average f value
            for (int i = 0; i < sys_size; i++)
            {
                ffavg_[i] = 0.5 * (ff1_[i] + ff2_[i]);
            }

            // Corrected value
            for (int i = 0; i < sys_size; i++)
            {
                yy2_[i] = yy1_[i] + (in_h * ffavg_[i]);
            }

            // Error calculation
            err = fabs((yy2_[0] - yy2pred_[0]) / (yy2_[0] + DBL_EPSILON));
            for (int i = 1; i < sys_size; ++i)
            {
                real temp_err = fabs((yy2_[i] - yy2pred_[i]) / (yy2_[i] + DBL_EPSILON));
                if (temp_err > err)
                {
                    err = temp_err;
                }
            }

            iter++;
        }

        if (iter > max_iter)
        {
            printf("\n Maximum number of iterations: %d\n", max_iter);
            printf("Exceeded at time: %f\n", t1);
            printf("In function '_dynorb_heun_'\n\n");
            break;
        }

        // Update
        t += in_h;
        memcpy(yy, yy2_, sys_size * sizeof(real));

        // if (step >= in_n_steps)
        // {
        //     out_tt = (real *)realloc(out_tt, (in_n_steps + 1) * sizeof(real));
        //     out_yyt = (real *)realloc(out_yyt, (in_n_steps + 1) * sys_size * sizeof(real));
        //     out_tt[step] = t;
        //     memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        //     printf("Number of estimated steps (in_n_steps = %d) surpassed. Using 'realloc'.\n", in_n_steps);
        // }
        // else
        // {
        //     out_tt[step] = t;
        //     memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        // }
        out_tt[step] = t;
        memcpy(&out_yyt[step * sys_size], yy, sys_size * sizeof(real));
        step++;
    }

    // Free Memory:
    free(ffavg_);
    free(yy2pred_);
    free(ff2_);
    free(ff1_);
    free(yy2_);
    free(yy1_);
    free(yy);
}

#endif // DYNORB_IMPLEMENTATION
