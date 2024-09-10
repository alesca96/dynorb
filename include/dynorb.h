/* ****************************************************
 * ****************** DYNORB HEADER *******************
 * **************************************************** */

/* ASSUMPTIONS:
 * 1. Code uses COLUMN MAJOR convention (similar to CBLAS).
 * 2. Strides other than 1 are not supported in this implementation.
 */

#ifndef DYNORB_H_
#define DYNORB_H_

/* SATNDARD C-LIBS: */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* CBLAS SWITCH: */
#ifdef USE_CBLAS
#include <cblas.h>
#endif

/* DOUBLE-FLOAT SWITCH: */
#ifdef USE_DOUBLE
typedef double real;
const bool real_is_double = 1;
#elif defined(USE_FLOAT)
typedef float real;
const bool real_is_double = 0;
#endif

/* MACROS: */
#define EL(M, m, n, i, j) M[(j * (m)) + i]

/* ODE FUNCTION AND SYSTEM: */
typedef void(_dynorb_odeFun)(const real t, const real *yy, const void *params, real *dyydt);
typedef struct
{
    _dynorb_odeFun *odeFunction; // Pointer to the ODE function
    const void *params;          // Pointer to the parameters for the ODE function
    const real *yy0;             // Pointer to the initial state array
    const real t0;               // Initial time
    const real t1;               // Final time
    const int sys_size;          // Size of the system (number of equations)
    real *tt;                    // Time Steps of the Solution
    real *YY_t;                  // Solution Array
} _dynorb_odeSys;

/* UTILITY FUNCTIONS: */
int _dynorb_min(int a, int b);
int _dynorb_max(int a, int b);
real _dynorb_rmin(real a, real b);
real _dynorb_rmax(real a, real b);
void _dynorb_rmprint(const real *A, const int m, const int n);

/* Basic Linear Algebra Subroutines: */
// Thin wrappers of CBLAS if USE_CBLAS. Else custum implementation.
void _dynorb_raxpy(const int n, const real alpha, const real *xx, real *yy);
void _dynorb_rgemv(const bool TransposeA, const int m, const int n, const real alpha, const real *A, const real *xx, const real beta, real *yy);
void _dynorb_rvcopy(const int n, const real *src_xx, real *dst_yy);
void _dynorb_rmcopy(const real *src_A, const int src_m, const int src_n, const int src_start_i, const int src_start_j, const int cpy_m, const int cpy_n, real *dst_B, const int dst_m, const int dst_n, const int dst_start_i, const int dst_start_j);

/* CORE FUNCTIONS: */

/**
 *
 * @brief Performs Euler (RK1) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Euler (Runge-Kutta order 1) method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] h Time step size for the integration.
 * @param[in] n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk1(_dynorb_odeSys *sys, const real h, const int n_steps);

/**
 *
 * @brief Performs Heun (RK2) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Heun (Runge-Kutta order 2) method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] h Time step size for the integration.
 * @param[in] n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk2(_dynorb_odeSys *sys, const real h, const int n_steps);

/**
 *
 * @brief Performs Heun (RK3) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Runge-Kutta order 3 method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] h Time step size for the integration.
 * @param[in] n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk3(_dynorb_odeSys *sys, const real h, const int n_steps);

/**
 *
 * @brief Performs Heun (RK4) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Runge-Kutta order 4 method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] h Time step size for the integration.
 * @param[in] n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk4(_dynorb_odeSys *sys, const real h, const int n_steps);

/**
 *
 * @brief Performs Heun Predictor-Corrector numerical integration for ODE systems.
 *
 * This function integrates an ODE system using Heun Predictor-Corrector method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] h Time step size for the integration.
 * @param[in] n_steps Number of steps for the integration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_heun_(_dynorb_odeSys *sys, real h, const int n_steps);

#endif // DYNORB_H_

/*
 * **************************************************** *
 * ************* DYNORN IMPLEMENTATION ************* *
 * **************************************************** *
 */

#ifdef DYNORB_IMPLEMENTATION

#ifdef USE_CBLAS /* CBLAS wrappers: */

void _dynorb_rvcopy(const int n, const real *src_xx, real *dst_yy)
{ // Vector copy: src_x->dst_y
    if (real_is_double)
    {
        cblas_dcopy(n, (const double *)src_xx, 1, (double *)dst_yy, 1);
    }
    else // Real is Float
    {
        cblas_scopy(n, (const float *)src_xx, 1, (float *)dst_yy, 1);
    }
}

void _dynorb_rmcopy(const real *src_A, const int src_M, const int src_N, const int src_StartRow, const int src_StartCol, const int cpyM, const int cpyN, real *dst_B, const int dst_M, const int dst_N, const int dst_StartRow, const int dst_StartCol)
{ // Matrix copy: A->B (also subportions), assumes CblasColMajor
    if (src_StartRow + cpyM > src_M || src_StartCol + cpyN > src_N)
    {
        printf("Error: Source submatrix dimensions exceed bounds of source matrix.\n");
        return;
    }
    if (dst_StartRow + cpyM > dst_M || dst_StartCol + cpyN > dst_N)
    {
        printf("Error: Destination submatrix dimensions exceed bounds of destination matrix.\n");
        return;
    }
    // Loop through each column in the source submatrix and copy it to the destination submatrix
    for (int col = 0; col < cpyN; ++col)
    {
        // Compute the offset to the first element of the column in both the source and destination
        const real *src_col = src_A + (src_StartCol + col) * src_M + src_StartRow;
        real *dst_col = dst_B + (dst_StartCol + col) * dst_M + dst_StartRow;

        // Copy the column (which is a vector) from source to destination
        _dynorb_rvcopy(cpyM, src_col, dst_col);
    }
}

void _dynorb_raxpy(const int n, const real alpha, const real *xx, real *yy)
{ // Vector Sum: yy = alpha*xx + yy

    if (real_is_double)
    {
        cblas_daxpy(n, (const double)alpha, (const double *)xx, 1, (double *)yy, 1);
    }
    else // Real is Float
    {
        cblas_saxpy(n, (const float)alpha, (const float *)xx, 1, (float *)yy, 1);
    }
}

void _dynorb_rgemv(const bool TransposeA, const int m, const int n, const real alpha, const real *A, const real *xx, const real beta, real *yy)
{ // matrix Vector Multiplication: yy = alpha*(A*xx) + beta*yy

    if (real_is_double)
    {
        if (TransposeA)
        {
            cblas_dgemv(CblasColMajor, CblasTrans, m, n, (const double)alpha, (const double *)A, m, (const double *)xx, 1, (const double)beta, (double *)yy, 1);
        }
        else // Not Transposed
        {
            cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, (const double)alpha, (const double *)A, m, (const double *)xx, 1, (const double)beta, (double *)yy, 1);
        }
    }
    else // Real is Float
    {
        if (TransposeA)
        {
            cblas_sgemv(CblasColMajor, CblasTrans, m, n, (const float)alpha, (const float *)A, m, (const float *)xx, 1, (const float)beta, (float *)yy, 1);
        }
        else // Not Transposed
        {
            cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, (const float)alpha, (const float *)A, m, (const float *)xx, 1, (const float)beta, (float *)yy, 1);
        }
    }
}

#endif

#ifndef USE_CBLAS /* BLAS-like Custom Implementation wrappers: */
/* TODO: custom implementation-> */
void IMPLEMENT_ALTERNATIVE_TO_CBLAS_FUNCTIONS(void)
{
    printf("I HAVE TO IMPLEMENT ALTERNATIVE TO CBLAS FUNCTIONS.\n");
}
#endif

/* Utility functions: */
int _dynorb_min(int a, int b)
{
    return (a < b) ? a : b;
}

int _dynorb_max(int a, int b)
{
    return (a > b) ? a : b;
}

real _dynorb_rmin(real a, real b)
{
    return (a < b) ? a : b;
}

real _dynorb_rmax(real a, real b)
{
    return (a > b) ? a : b;
}

void _dynorb_rmprint(const real *A, const int m, const int n)
{
    printf("[\n");
    for (int i = 0; i < m; ++i)
    {
        printf("    ");
        for (int j = 0; j < n; ++j)
        {
            printf("%.4f", EL(A, m, n, i, j));
            if (j < n - 1)
            {
                printf(", ");
            }
        }
        printf(";\n");
    }
    printf("]\n");
}

/* Core Library Functions: */
void _dynorb_rk1(_dynorb_odeSys *sys, const real h, const int n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = sys->params;               // Pointer to params
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

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
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives in [t, t+h]
        real t_inner = t + a * h;
        odeFunction(t_inner, yy, params, ff);

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += h * c * ff[k];
            }
        }

        // Update Integration Time:
        t += h;

        // Store results:
        if (t > t1)
        {
            memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
            tt[step] = t;
            printf("_dynorb_rk1 : Breaking From Loop<t = %f [s]>", t);
        }
        memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
        tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk1: Begin Freeing Memory:\n");
    free(ff);
    free(yy);
    // printf("_dynorb_rk1: Done Freeing Memory:\n");
}

void _dynorb_rk2(_dynorb_odeSys *sys, const real h, const int n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = sys->params;               // Pointer to params
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

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
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * h;
            memcpy(yy_inner, yy, sys_size * sizeof(real));
            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, params, &ff[i * sys_size]);
        }

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += h * c[i] * ff[i * sys_size + k];
            }
        }

        // Update Integration Time:
        t += h;

        // Store results:
        if (t > t1)
        {
            memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
            tt[step] = t;
            printf("_dynorb_rk2 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
        tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk2: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    // printf("_dynorb_rk2: Done Freeing Memory:\n");
}

void _dynorb_rk3(_dynorb_odeSys *sys, const real h, const int n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = sys->params;               // Pointer to params
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

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
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * h;
            memcpy(yy_inner, yy, sys_size * sizeof(real));
            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, params, &ff[i * sys_size]);
        }

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += h * c[i] * ff[i * sys_size + k];
            }
        }

        // Update Integration Time:
        t += h;

        // Store results:
        if (t > t1)
        {
            memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
            tt[step] = t;
            printf("_dynorb_rk3 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
        tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk3: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    // printf("_dynorb_rk3: Done Freeing Memory:\n");
}

void _dynorb_rk4(_dynorb_odeSys *sys, const real h, const int n_steps)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = sys->params;               // Pointer to params
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

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
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * h;
            memcpy(yy_inner, yy, sys_size * sizeof(real));
            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, params, &ff[i * sys_size]);
        }

        // Update the state:
        for (int k = 0; k < sys_size; ++k)
        {
            for (int i = 0; i < n_stages; ++i)
            {
                yy[k] += h * c[i] * ff[i * sys_size + k];
            }
        }

        // Update Integration Time:
        t += h;

        // Store results:
        if (t > t1)
        {
            memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
            tt[step] = t;
            printf("_dynorb_rk4 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
        tt[step] = t;
    }

    // Free Memory:
    // printf("_dynorb_rk4: Begin Freeing Memory:\n");
    free(yy_inner);
    free(ff);
    free(yy);
    // printf("_dynorb_rk4: Done Freeing Memory:\n");
}

void _dynorb_heun_(_dynorb_odeSys *sys, real h, const int n_steps)
{
    // Open up _dynorb_odeSys
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *params = sys->params;               // Pointer to params
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Tolerance and _dynorb_Max number of steps:
    const real tol = 1.e-6;
    const int rmax_iter = 100;

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
    tt[0] = t0;
    memcpy(YY_t, yy0, sys_size * sizeof(real));

    // Allocate Memory for state and derivatives at interval boundaries:
    real *yy1_ = (real *)malloc(sys_size * sizeof(real));
    real *yy2_ = (real *)malloc(sys_size * sizeof(real));
    real *ff1_ = (real *)malloc(sys_size * sizeof(real));
    real *ff2_ = (real *)malloc(sys_size * sizeof(real));
    real *yy2pred_ = (real *)malloc(sys_size * sizeof(real));
    real *ffavg_ = (real *)malloc(sys_size * sizeof(real));

    // Main Loop
    int step = 1;
    while (t < t1 && step < n_steps)
    {
        // Step Size:
        h = _dynorb_rmin(h, t1 - t);

        // Left Boundary of Interval:
        real t1_ = t;
        memcpy(yy1_, yy, sys_size * sizeof(real));
        odeFunction(t1_, yy1_, params, ff1_);

        // Compute yy2_
        for (int i = 0; i < sys_size; ++i)
        {
            yy2_[i] = yy1_[i] + ff1_[i] * h;
        }

        // Right Boundary of the Interval:
        real t2_ = t1_ + h;
        real err = tol + 1;
        int iter = 0;

        // Predictor-Corrector Loop
        while (err > tol && iter <= rmax_iter)
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
                yy2_[i] = yy1_[i] + (h * ffavg_[i]);
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

        if (iter > rmax_iter)
        {
            printf("\n Maximum number of iterations: %d\n", rmax_iter);
            printf("Exceeded at time: %f\n", t1);
            printf("In function '_dynorb_heun_'\n\n");
            break;
        }

        // Update
        t += h;
        memcpy(yy, yy2_, sys_size * sizeof(real));
        /*
        if (step >= n_steps)
         {
             tt = (real *)realloc(tt, (n_steps + 1) * sizeof(real));
             YY_t = (real *)realloc(YY_t, (n_steps + 1) * sys_size * sizeof(real));
             tt[step] = t;
             memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
             printf("Number of estimated steps (n_steps = %d) surpassed. Using 'realloc'.\n", n_steps);
         }
         else
         {
             tt[step] = t;
             memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
         }
        */
        tt[step] = t;
        memcpy(&YY_t[step * sys_size], yy, sys_size * sizeof(real));
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
