/* ****************************************************
 * ****************** DYNORB HEADER *******************
 * **************************************************** */

/* ASSUMPTIONS:
 * 1. Code uses COLUMN MAJOR convention (see CBLAS).
 * 2. Strides other than 1 are not supported in this implementation.
 * 3. Type 'real' can only be float (single precision) or double (double precision).
 */

/* TODO: Short/Medium Term
 * 1. Implement Curtis Algorithms using wrappers of CBLAS and LAPACKE
 * 2. Reduce (completly eliminate when possible) heap usage
 * 3. Respect as much as possible NASA POWER OF TEN
 * 4. Implement all new functions as: _dynorb_rcamelCaseCaseCase (r-->real)
 */

/* TODO: Long Term
 * 1. Implement Custom subroutines to substitute CBLAS and LAPACKE
 * 2. Completly eliminate heap usage (if impossible, remove the function from library)
 * 3. Enforce NASA POWER OF TEN
 * 4. Conver all functions as: _dynorb_rcamelCaseCaseCase (r-->real)
 */

/* NASA POWER OF TEN RULES (http://web.eecs.umich.edu/~imarkov/10rules.pdf):
 * 1. Avoid complex flow constructs, such as goto and recursion.
 * 2. All loops must have fixed bounds. This prevents runaway code.
 * 3. Avoid heap memory allocation.
 * 4. Restrict functions to a single printed page.
 * 5. Use a minimum of two runtime assertions per function.
 * 6. Restrict the scope of data to the smallest possible.
 * 7. Check the return value of all non-void functions, or cast to void to indicate the return value is useless.
 * 8. Use the preprocessor sparingly.
 * 9. Limit pointer use to a single dereference, and do not use function pointers.
 * 10. Compile with all possible warnings active; all warnings should then be addressed before release of the software.
 */

#ifndef DYNORB_H_
#define DYNORB_H_

/* ====================================================================================================
 * SATNDARD C-LIBS:
 * ==================================================================================================== */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ====================================================================================================
 * CBLAS SWITCH:
 * ==================================================================================================== */
#ifdef USE_CBLAS
#include <cblas.h>
#endif

/* ====================================================================================================
 * DOUBLE-FLOAT SWITCH:
 * ==================================================================================================== */
#ifdef USE_DOUBLE
typedef double real;
const bool real_is_double = 1;
#elif defined(USE_FLOAT)
typedef float real;
const bool real_is_double = 0;
#endif

/* ====================================================================================================
 * DYNORB MACROS:
 * ==================================================================================================== */
#define _dynorb_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#define _dynorb_R_E 6378.1370     // [km]
#define _dynorb_MU_E 398600.44188 // Earth grav. parameter [km^3 /s^2]

/* ====================================================================================================
 * DYNORB GLOBAL VARIABLES:
 * ==================================================================================================== */

/* ====================================================================================================
 * DYNORB CORE FYNCTION TYPES:
 * ==================================================================================================== */
typedef real (*_dynorb_nonLinScalFun)(real a, void *funParams);                                 // Non-linear function
typedef void(_dynorb_odeFun)(const real t, const real *yy, const void *odeParams, real *dyydt); // ODE-function

/* ====================================================================================================
 * DYNORB CORE STRUCTURES:
 * ==================================================================================================== */

typedef struct
{
    _dynorb_odeFun *odeFunction; // Pointer to the ODE function
    void *odeParams;             // Pointer to the parameters (structure) for the ODE function
    int sys_size;                // Size of the system (number of equations)
    real t0;                     // Initial time
    real tf;                     // Final time
    real *yy0;                   // Pointer to the initial state array
    real *tt;                    // Time Steps of the Solution
    real *YY_t;                  // Solution Array
} _dynorb_odeSys;

typedef struct
{
    real t0;     // Initial time
    real tf;     // Final time
    real h;      // Initial Step Size
    int n_steps; // Initial Number of steps
} _dynorb_solverConf;

/* ====================================================================================================
 * DYNORB USER SUPPORT STRUCTURES:
 * ==================================================================================================== */

typedef struct
{
    real m1; // Body 1 mass [kg]
    real m2; // Body 2 mass [kg]
    real G;  // Universal gravit. constant [km^3/kg/s^2]

} _dynorb_twoBodyAbsParams;

typedef struct
{
    real m;  // Body Mass [kg]
    real mu; // Body grav. parameter [km^3 /s^2]

} _dynorb_twoBodyRelParams;

typedef struct
{
    real mu1;  // Body 1 grav. parameter [km^3 /s^2]
    real mu2;  // Body 2 grav. parameter [km^3 /s^2]
    real r12;  // Body 1 - Body 2 distance [km]
    real pi_1; // Body 1 mass ratio [-]: pi_1 = m1/(m1+m2)
    real pi_2; // Body 2 mass ratio [-]: pi_2 = m2/(m1+m2)
    real W;    // Angular velocity Body 2 around Body 1 [rad/s]:  W = sqrt((mu1+mu2)) / (r12 * r12 * r12));
    real x1;   // x-coordinates of the Body 1 relative to the Body 1 - Body 2 barycenter [km]: x1 = -1.0 * pi_2 * r12;
    real x2;   // x-coordinates of the Body 2 relative to the Body 1 - Body 2 barycenter [km]: x2 = 1.0 * pi_1 * r12;

} _dynorb_threeBodyRestrictParams;

/* ==========================================================================================================
 * DYNORB BASIC LINEAR ALGEBRA SUBROUTINES:
 * ========================================================================================================== */

real _dynorb_rdot(const int n, const real *xx, const real *yy);
real _dynorb_rnrm2(const int n, const real *xx);
void _dynorb_rscal(const int n, const real alpha, real *xx);
void _dynorb_raxpy(const int n, const real alpha, const real *xx, real *yy);
void _dynorb_rgemv(const bool TransposeA, const int m, const int n, const real alpha, const real *A, const real *xx, const real beta, real *yy);
void _dynorb_rvcopy(const int n, const real *src_xx, real *dst_yy);
void _dynorb_rmcopy(const real *src_A, const int src_m, const int src_n, const int src_start_i, const int src_start_j, const int cpy_m, const int cpy_n, real *dst_B, const int dst_m, const int dst_n, const int dst_start_i, const int dst_start_j);

/* ====================================================================================================
 * DYNORB UTILITIES FUNCTIONS:
 * ==================================================================================================== */

real _dynorb_rad2deg(real angle_rad);
real _dynorb_deg2rad(real angle_deg);
real _dynorb_max_abs_component(int n, const real *xx);
real _dynorb_eps(real a);
int _dynorb_min(int a, int b);
int _dynorb_max(int a, int b);
real _dynorb_rmin(real a, real b);
real _dynorb_rmax(real a, real b);
void _dynorb_rmprint(const real *A, const int m, const int n);
void _dynorb_rvfill(real *xx, const int n, const real a);
static inline real _dynorb_rel(const real *A, int m, int i, int j);
void _dynorb_rcol(const real *A, int m, int j, real *column);
void _dynorb_rrow(const real *A, int m, int n, int i, real *row);

/* ====================================================================================================
 * DYNORB USER SUPPORT FUNCTIONS:
 * ==================================================================================================== */

void _dynorb_free(_dynorb_odeSys *ode_system);
void _dynorb_configure_dynamic(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration, _dynorb_odeFun *ode_function, void *odeParams, real *yy0, const int sys_size, const real t0, const real tf, const real h);
void _dynorb_configure_static(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration, _dynorb_odeFun *ode_function, void *odeParams, real *yy0, const int sys_size, const real t0, const real tf, const real h);

/* ====================================================================================================
 * DYNORB NUMERICS:
 * ==================================================================================================== */

real _dynorb_rstumpffS(const real z);
real _dynorb_rstumpffC(const real z);
real _dynorb_rbisect(_dynorb_nonLinScalFun function, void *funParams, real a, real b, const real tol);
void _dynorb_rrk1(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);
void _dynorb_rrk2(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);
void _dynorb_rrk3(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);
void _dynorb_rrk4(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);
void _dynorb_rheun(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration, const real tol, const int max_iter);

/* ====================================================================================================
 * DYNORB ORBITAL MECHANICS FUNCTIONS:
 * ==================================================================================================== */

void _dynorb_twoBodyAbsFun(const real t, const real *yy, const void *params, real *ff);
void _dynorb_twoBodyRelFun(const real t, const real *yy, const void *params, real *ff);
void _dynorb_threeBodyRestrictFun(const real t, const real *yy, const void *params, real *ff);
real _dynorb_rkeplerH(real e, real M_hyperbola, real tol);
real _dynorb_rkeplerE(real e, real M_ellipse, real tol);
void _dynorb_rlagrangeCoefficientsFrom_yy0_Dth(const real mu, const real Dth, const real *yy0, real *fgdfdg);
void _dynorb_yy_From_yy0_Dth(const real mu, const real Dth, const real *yy0, real *yy);
real _dynorb_rkeplerUniversal(const real mu, const real Dt, const real r0, const real vr0, const real alpha, const real tol, const int max_iter);
void _dynorb_rfgFromChi(const real mu, const real Dt, const real chi, const real r0, const real alpha, real *fgdfdg);
void _dynorb_rdfdgFromChi(const real mu, const real chi, const real r, const real r0, const real alpha, real *fgdfdg);
void _dynorb_ryy_From_yy0_Dt(const real mu, const real Dt, const real *yy0, real *yy, const real tol, const int max_iter);

#endif // DYNORB_H_
/*
 *
 *
 *
 *
 *
 *
 *
 */
#ifdef DYNORB_IMPLEMENTATION

/* ==========================================================================================================
 * DYNORB BASIC LINEAR ALGEBRA SUBROUTINES:
 * ========================================================================================================== */
#ifdef USE_CBLAS // --> CBLAS wrappers

void _dynorb_rvcopy(const int n, const real *src_xx, real *dst_yy)
{ // LEVEL 1: Vector copy: src_x->dst_y
    if (real_is_double)
    {
        cblas_dcopy(n, (const double *)src_xx, 1, (double *)dst_yy, 1);
    }
    else // Real is Float
    {
        cblas_scopy(n, (const float *)src_xx, 1, (float *)dst_yy, 1);
    }
}

void _dynorb_rscal(const int n, const real alpha, real *xx)
{ // LEVEL 1: Scale a vectory by constant: xx = alpha*xx
    if (real_is_double)
    {
        cblas_dscal(n, (const double)alpha, (double *)xx, 1);
    }
    else // Real is Float
    {
        cblas_sscal(n, (const float)alpha, (float *)xx, 1);
    }
}

void _dynorb_raxpy(const int n, const real alpha, const real *xx, real *yy)
{ // LEVEL 1: Vector Sum: yy = alpha*xx + yy

    if (real_is_double)
    {
        cblas_daxpy(n, (const double)alpha, (const double *)xx, 1, (double *)yy, 1);
    }
    else // Real is Float
    {
        cblas_saxpy(n, (const float)alpha, (const float *)xx, 1, (float *)yy, 1);
    }
}

real _dynorb_rdot(const int n, const real *xx, const real *yy)
{ // LEVEL 1: Dot Product: dot =  xx.T*yy. Returns a (real) scalar value, inputs not modified.
    if (real_is_double)
    {
        double dot = cblas_ddot(n, (double *)xx, 1, (double *)yy, 1);
        return dot;
    }
    else // Real is Float
    {
        float dot = cblas_sdot(n, (float *)xx, 1, (float *)yy, 1);
        return dot;
    }
}

real _dynorb_rnrm2(const int n, const real *xx)
{ // LEVEL 1: Norm 2: nrm2 =  xx.T*xx. Returns a (real) scalar value, input not modified.
    if (real_is_double)
    {
        double nrm2 = cblas_dnrm2(n, (double *)xx, 1);
        return nrm2;
    }
    else // Real is Float
    {
        float nrm2 = cblas_snrm2(n, (float *)xx, 1);
        return nrm2;
    }
}

void _dynorb_rmcopy(const real *src_A, const int src_M, const int src_N,
                    const int src_StartRow, const int src_StartCol,
                    const int cpy_M, const int cpy_N,
                    real *dst_B, const int dst_M, const int dst_N,
                    const int dst_StartRow, const int dst_StartCol)
{ // LEVEL 2: Matrix copy: A->B (also subportions), assumes CblasColMajor
    if (src_StartRow + cpy_M > src_M || src_StartCol + cpy_N > src_N)
    {
        printf("Error: Source submatrix dimensions exceed bounds of source matrix.\n");
        return;
    }
    if (dst_StartRow + cpy_M > dst_M || dst_StartCol + cpy_N > dst_N)
    {
        printf("Error: Destination submatrix dimensions exceed bounds of destination matrix.\n");
        return;
    }
    // Loop through each column in the source submatrix and copy it to the destination submatrix
    for (int col = 0; col < cpy_N; ++col)
    {
        // Compute the offset to the first element of the column in both the source and destination
        const real *src_col = src_A + (src_StartCol + col) * src_M + src_StartRow;
        real *dst_col = dst_B + (dst_StartCol + col) * dst_M + dst_StartRow;

        // Copy the column (which is a vector) from source to destination
        _dynorb_rvcopy(cpy_M, src_col, dst_col);
    }
}

void _dynorb_rgemv(const bool TransposeA, const int m, const int n, const real alpha, const real *A, const real *xx, const real beta, real *yy)
{ // LEVEL 2: Matrix Vector Multiplication: yy = alpha*(A*xx) + beta*yy

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

#else // #ifndef USE_CBLAS // --> BLAS-like TODO: Custom Implementation
void IMPLEMENT_ALTERNATIVE_TO_CBLAS_FUNCTIONS(void)
{
    printf("I HAVE TO IMPLEMENT ALTERNATIVE TO CBLAS FUNCTIONS.\n");
}
#endif

/* ====================================================================================================
 * DYNORB UTILITIES FUNCTIONS:
 * ==================================================================================================== */
real _dynorb_rad2deg(real angle_rad)
{
    real angle_deg = (angle_rad / _dynorb_PI) * 180.0;
    return angle_deg;
}

real _dynorb_deg2rad(real angle_deg)
{
    real angle_rad = (angle_deg / 180.0) * _dynorb_PI;
    return angle_rad;
}

real _dynorb_max_abs_component(int n, const real *xx)
{ // Function to compute the maximum absolute component of a xxtor
    real max_abs = fabs(xx[0]);
    for (int i = 1; i < n; ++i)
    {
        real tmp = (real)fabs(xx[i]);
        if (tmp > max_abs)
        {
            max_abs = tmp;
        }
    }

    return max_abs;
}

real _dynorb_eps(real a)
{
    real eps = 1.0;

    // Iteratively divide eps by 2 until x + eps equals x
    while (a + eps != a)
    {
        eps /= 2.0;
    }

    // Return the previous value of eps (smallest such that x + eps != x)
    return eps * 2.0;
}

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
            printf("%.5f", _dynorb_rel(A, m, i, j));
            if (j < n - 1)
            {
                printf(", ");
            }
        }
        printf(";\n");
    }
    printf("]\n");
}

void _dynorb_rvfill(real *xx, const int n, const real a)
{ // Fills vector components xx[i] = a for i  = [0, n)
    for (int i = 0; i < n; ++i)
    {
        xx[i] = a;
    }
}

static inline real _dynorb_rel(const real *A, int m, int i, int j)
{
    return A[(j * m) + i];
}

void _dynorb_rcol(const real *A, int m, int j, real *column)
{
    for (int i = 0; i < m; i++)
    {
        column[i] = _dynorb_rel(A, m, i, j); // Accessing element A[i, j]
    }
}

void _dynorb_rrow(const real *A, int m, int n, int i, real *row)
{
    for (int j = 0; j < n; j++)
    {
        row[j] = _dynorb_rel(A, m, i, j); // Accessing element A[i, j]
    }
}

/* ====================================================================================================
 * DYNORB USER SUPPORT FUNCTIONS:
 * ==================================================================================================== */

void _dynorb_free(_dynorb_odeSys *ode_system)
{
    free(ode_system->tt);
    free(ode_system->YY_t);
}

void _dynorb_configure_dynamic(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration, _dynorb_odeFun *ode_function, void *odeParams, real *yy0, const int sys_size, const real t0, const real tf, const real h)
{
    // Compute (initial) number of steps:
    int n_steps = ceil(((tf - t0) / h)); // int n_steps = (int)((((tf - t0)) / h) + 1.0);

    // Allocate memory for the time steps and solution array
    ode_system->tt = (real *)calloc(n_steps, sizeof(real));
    ode_system->YY_t = (real *)calloc(n_steps, sys_size * sizeof(real));

    // Check for successful allocation
    if (ode_system->tt == NULL || ode_system->YY_t == NULL)
    {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        if (ode_system->tt != NULL)
            free(ode_system->tt); // Free memory if tt was allocated
        exit(EXIT_FAILURE);
    }

    // Print dynamic memory allocation information
    printf("\nDynamic Memory Allocation:\n");
    printf("Remember to free the memory after usage:\n");
    printf("    _dynorb_free(&ode_system).\n\n");

    // Print stack memory allocation instructions
    printf("For Stack Memory Allocation, use _dynorb_configure_static() and add these lines to your code:\n");
    printf("    real tt[solver_configuration.n_steps];\n");
    printf("    real YY_t[solver_configuration.n_steps * ode_system.sys_size];\n");
    printf("    ode_system.tt = tt;\n");
    printf("    ode_system.YY_t = YY_t;\n\n");

    // Copy values in ode_system structure:
    ode_system->odeFunction = ode_function;
    ode_system->odeParams = odeParams;
    ode_system->sys_size = sys_size;
    ode_system->t0 = t0;
    ode_system->tf = tf;
    ode_system->yy0 = yy0;

    // Copy values in solver_configuration structure:
    solver_configuration->t0 = t0;
    solver_configuration->tf = tf;
    solver_configuration->h = h;
    solver_configuration->n_steps = n_steps;

    // Print:
    printf("\nTime Step Size: <h = %f [s]>\n", solver_configuration->h);
    printf("Number of Steps: <solv_conf.n_steps = %d>\n", solver_configuration->n_steps);
}

void _dynorb_configure_static(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration, _dynorb_odeFun *ode_function, void *odeParams, real *yy0, const int sys_size, const real t0, const real tf, const real h)
{
    // Compute (initial) number of steps:
    int n_steps = ceil(((tf - t0) / h)); // int n_steps = (int)((((tf - t0)) / h) + 1.0);

    // Print stack memory allocation instructions
    ode_system->tt = NULL;
    ode_system->YY_t = NULL;

    // Copy values in ode_system structure:
    ode_system->odeFunction = ode_function;
    ode_system->odeParams = odeParams;
    ode_system->sys_size = sys_size;
    ode_system->t0 = t0;
    ode_system->tf = tf;
    ode_system->yy0 = yy0;

    // Copy values in solver_configuration structure:
    solver_configuration->t0 = t0;
    solver_configuration->tf = tf;
    solver_configuration->h = h;
    solver_configuration->n_steps = n_steps;

    // Print:
    printf("\nTime Step Size: <h = %f [s]>\n", solver_configuration->h);
    printf("Number of Steps: <solv_conf.n_steps = %d>\n", solver_configuration->n_steps);
}

/* ====================================================================================================
 * DYNORB NUMERICS:
 * ==================================================================================================== */

real _dynorb_rstumpffS(const real z)
{
    real S;
    if (z > 0)
    {
        real sqrt_z = sqrt(z);
        S = (sqrt_z - sin(sqrt_z)) / pow(sqrt_z, 3.0);
    }
    else if (z < 0)
    {
        real sqrt_z = sqrt(-z);
        S = (sinh(sqrt_z) - sqrt_z) / pow(sqrt_z, 3.0);
    }
    else
    {
        S = 1.0 / 6.0;
    }
    return S;
}

real _dynorb_rstumpffC(const real z)
{
    real C;
    if (z > 0)
    {
        C = (1 - cos(sqrt(z))) / z;
    }
    else if (z < 0)
    {
        C = (cosh(sqrt(-z)) - 1) / (-z);
    }
    else
    {
        C = 0.5;
    }
    return C;
}

real _dynorb_rbisect(_dynorb_nonLinScalFun function, void *funParams, real a, real b, const real tol)
{ // Implementation of Bisection Method
    // Midpoit
    real c = 0;
    // Number of Iteration: TODO: fix fabs, works inly with double
    real n = ceil(log(fabs(b - a) / tol) / log(2));
    for (int i = 0; i < n; ++i)
    {
        c = 0.5 * (a + b); // midpoint update
        real f_a = function(a, funParams);
        real f_c = function(c, funParams);
        if (f_a * f_c > 0)
        {
            a = c;
        }
        else
        {
            b = c;
        }
    }
    return c;
}

void _dynorb_rrk1(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real tf = sys->tf;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Open up _dynorb_solvConf [TODO: remove this]
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    //  RK1 (Euler)
    const int n_stages = 1;

    // Allocate Memory for Current State:
    real yy[sys_size];
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real ff[sys_size * n_stages]; // (ff = dyy/dt)

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at time instant t:
        odeFunction(t, yy, odeParams, ff);

        // Update state: Vector Sum: yy = alpha*xx + yy
        _dynorb_raxpy(sys_size, h, ff, yy);

        // Update Integration Time:
        t = (step + 1) * h; // t += h;

        // Store results:
        _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
        tt[step] = t;

        if (t > tf)
        {
            printf("\n_dynorb_rrk1: Breaking From Loop at <t = %f [s]>", t);
            break;
        }
    }
}

void _dynorb_rrk2(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up structures [TODO: remove these copies]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real tf = sys->tf;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size];
    real yy_inner[sys_size];
    real ff1_[sys_size];
    real ff2_[sys_size];
    real ff1_plus_ff2_[sys_size];

    // Integration Time Instant:
    real t = t0;
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at time 1 = t
        odeFunction(t, yy, odeParams, ff1_);
        // Evaluate Time Derivatives at time 2 = t+h
        _dynorb_rvcopy(sys_size, yy, yy_inner);
        _dynorb_raxpy(sys_size, h, ff1_, yy_inner);
        odeFunction((t + h), yy_inner, odeParams, ff2_);
        // Update state yy:
        _dynorb_rscal(sys_size, 0.5, ff1_);
        _dynorb_rvcopy(sys_size, ff1_, ff1_plus_ff2_);
        _dynorb_raxpy(sys_size, 0.5, ff2_, ff1_plus_ff2_);
        _dynorb_raxpy(sys_size, h, ff1_plus_ff2_, yy);
        // Update Integration Time:
        t = (step + 1) * h; // t += h;
        // Store results:
        _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
        tt[step] = t;
        // Safety check:
        if (t > tf)
        {
            printf("\n_dynorb_rrk2: Breaking From Loop at <t = %f [s]>", t);
            break;
        }
    }
}

void _dynorb_rrk3(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up structures [TODO: remove these copies]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real tf = sys->tf;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size];
    real yy_inner[sys_size];
    real ff1_[sys_size];
    real ff2_[sys_size];
    real ff3_[sys_size];
    real two_ff2_[sys_size];
    real minus_ff1_plus_two_ff2_[sys_size];

    // Integration Time Instant:
    real t = t0;
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at time 1 = t
        odeFunction(t, yy, odeParams, ff1_);
        // Evaluate Time Derivatives at time 2 = t+(0.5*h)
        _dynorb_rvcopy(sys_size, yy, yy_inner);
        _dynorb_raxpy(sys_size, (0.5 * h), ff1_, yy_inner);
        odeFunction((t + (0.5 * h)), yy_inner, odeParams, ff2_);
        // Evaluate Time Derivatives at time 3 = t+h
        _dynorb_rvcopy(sys_size, yy, yy_inner);
        _dynorb_rvcopy(sys_size, ff2_, two_ff2_);
        _dynorb_rscal(sys_size, 2.0, two_ff2_);
        _dynorb_rvcopy(sys_size, two_ff2_, minus_ff1_plus_two_ff2_);
        _dynorb_raxpy(sys_size, -1.0, ff1_, minus_ff1_plus_two_ff2_);
        _dynorb_raxpy(sys_size, h, minus_ff1_plus_two_ff2_, yy_inner);
        odeFunction((t + h), yy_inner, odeParams, ff3_);
        // Update state yy:
        _dynorb_rscal(sys_size, (1.0 / 3.0), two_ff2_);
        _dynorb_raxpy(sys_size, (1.0 / 6.0), ff3_, two_ff2_);
        _dynorb_raxpy(sys_size, (1.0 / 6.0), ff1_, two_ff2_);
        _dynorb_raxpy(sys_size, h, two_ff2_, yy);

        // Update Integration Time:
        t = (step + 1) * h; // t += h;

        // Store results:
        _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
        tt[step] = t;

        if (t > tf)
        {
            printf("\n_dynorb_rrk3: Breaking From Loop at <t = %f [s]>", t);
            break;
        }
    }
}

void _dynorb_rrk4(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up structures [TODO: remove these copies]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real tf = sys->tf;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size];
    real yy_inner[sys_size];
    real ff1_[sys_size];
    real ff2_[sys_size];
    real ff3_[sys_size];
    real ff4_[sys_size];

    // Integration Time Instant:
    real t = t0;
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at time 1 = t
        odeFunction(t, yy, odeParams, ff1_);
        // Evaluate Time Derivatives at time 2 = t+(0.5*h)
        _dynorb_rvcopy(sys_size, yy, yy_inner);
        _dynorb_raxpy(sys_size, (0.5 * h), ff1_, yy_inner);
        odeFunction((t + (0.5 * h)), yy_inner, odeParams, ff2_);
        // Evaluate Time Derivatives at time 3 = t+(0.5*h)
        _dynorb_rvcopy(sys_size, yy, yy_inner);
        _dynorb_raxpy(sys_size, (0.5 * h), ff2_, yy_inner);
        odeFunction((t + (0.5 * h)), yy_inner, odeParams, ff3_);
        // Evaluate Time Derivatives at time 4 = t+h
        _dynorb_rvcopy(sys_size, yy, yy_inner);
        _dynorb_raxpy(sys_size, h, ff3_, yy_inner);
        odeFunction((t + h), yy_inner, odeParams, ff4_);
        // Update state yy:
        _dynorb_rscal(sys_size, (1.0 / 6.0), ff4_);
        _dynorb_raxpy(sys_size, (1.0 / 3.0), ff3_, ff4_);
        _dynorb_raxpy(sys_size, (1.0 / 3.0), ff2_, ff4_);
        _dynorb_raxpy(sys_size, (1.0 / 6.0), ff1_, ff4_);
        _dynorb_raxpy(sys_size, h, ff4_, yy);

        // Update Integration Time:
        t = (step + 1) * h; // t += h;

        // Store results:
        _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
        tt[step] = t;

        if (t > tf)
        {
            printf("\n_dynorb_rrk3: Breaking From Loop at <t = %f [s]>", t);
            break;
        }
    }
}

void _dynorb_rheun(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration, const real tol, const int max_iter)
{
    // Extract pointers from the input structs
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    const real t0 = sys->t0;                        // Initial time
    const real tf = sys->tf;                        // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array
    const real h = solver_configuration->h;
    int n_steps = solver_configuration->n_steps;

    // Internal predictor-corrector loop:
    real ee[sys_size]; // error vector
    real err = tol + 1.0;
    int iter = 0;

    // Allocate memory for current and predicted state:
    real t1_;
    real t2_;
    real yy[sys_size];
    real yy1_[sys_size];
    real yy2_[sys_size];
    real yy2p_[sys_size];
    real ff1_[sys_size];
    real ff2_[sys_size];
    real ffavg_[sys_size];

    // Set initial state:
    real t = t0;
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Numerical integration:
    for (int step = 0; step < n_steps; ++step)
    {
        t1_ = t;
        _dynorb_rvcopy(sys_size, yy, yy1_);
        odeFunction(t1_, yy1_, odeParams, ff1_);

        // Predictor step
        _dynorb_rvcopy(sys_size, yy1_, yy2_);
        _dynorb_raxpy(sys_size, h, ff1_, yy2_);
        t2_ = t1_ + h;

        err = tol + 1.0;
        iter = 0;

        // Corrector loop
        while (err > tol && iter <= max_iter)
        {
            _dynorb_rvcopy(sys_size, yy2_, yy2p_);    // Copy the predicted state
            odeFunction(t2_, yy2p_, odeParams, ff2_); // Compute f(t+h, y_pred)

            // Compute average of slopes
            _dynorb_rvcopy(sys_size, ff2_, ffavg_);
            _dynorb_raxpy(sys_size, 1.0, ff1_, ffavg_); // ffavg = ff1 + ff2
            _dynorb_rscal(sys_size, 0.5, ffavg_);       // ffavg = 0.5 * (ff1 + ff2)

            // Update the predicted state using the average slopes
            _dynorb_rvcopy(sys_size, yy1_, yy2_); // yy2_ = yy1_ + h * ffavg_
            _dynorb_raxpy(sys_size, h, ffavg_, yy2_);

            // Compute the error: ee = yy2_ - yy2p_
            _dynorb_rvcopy(sys_size, yy2_, ee); // ee = yy2_ - yy2p_
            _dynorb_raxpy(sys_size, -1.0, yy2p_, ee);

            // Compute the maximum relative error
            err = 0.0;
            for (int i = 0; i < sys_size; ++i)
            {
                real tmp = (real)fabs(ee[i] / (yy2_[i] + DBL_EPSILON));
                if (tmp > err)
                {
                    err = tmp;
                }
            }

            ++iter;
        }

        // Safety check for exceeding max iterations
        if (iter >= max_iter)
        {
            printf("\nMaximum number of iterations: %d\n", max_iter);
            printf("Exceeded at time: %f\n", t);
            printf("In function '_dynorb_rheun'. \n\n");
            break;
        }

        // Update integration time: t += h
        t += h;

        // Store the results:
        _dynorb_rvcopy(sys_size, yy2_, yy); // Update the state
        tt[step] = t;
        _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);

        // Stop if the final time is exceeded
        if (t >= tf)
        {
            printf("\n_dynorb_rheun: Breaking from loop at t = %f [s]\n", t);
            break;
        }
    }
}

void _dynorb_rrkf45(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration, const real tol)
{
    // Print WARNING:
    printf("\nWARNING: Arrays sys->tt and sys->YY_t must be dynamically allocated to use this function.\n");
    printf("         Please make sure to configure using _dynorb_configure_dynamic.\n");
    printf("         If this is the case already, ignore this warning.\n\n");

    // Open up structures [TODO: remove these copies]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real tf = sys->tf;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array
    real h = solver_configuration->h;               // Initial step
    int n_steps = solver_configuration->n_steps;

    // Number of stages:
    int n_stages = 6;

    // Coefficients:
    real a[6] = {0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0};
    real b[36] = {
        0.0, 1.0 / 4.0, 3.0 / 32.0, 1932.0 / 2197.0, 439.0 / 216.0, -8.0 / 27.0, // first column
        0.0, 0.0, 9.0 / 32.0, -7200.0 / 2197.0, -8.0, 2.0,                       // second column
        0.0, 0.0, 0.0, 7296.0 / 2197.0, 680.0 / 513.0, -3544.0 / 2565.0,         // third column
        0.0, 0.0, 0.0, 0.0, -845.0 / 4104.0, 1859.0 / 4104.0,                    // fourth column
        0.0, 0.0, 0.0, 0.0, 0.0, -11.0 / 40.0,                                   // fifth column
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0                                             // sixth column
    };
    real c_star[6] = {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0};
    real c[6] = {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size]; // High order solution (RK5)
    // real yy_star[sys_size]; // Low order solution (RK4)
    real yyi_[sys_size];
    real yy_inner[sys_size];
    real ee[sys_size]; // Truncation error vector
    real FF_[sys_size * n_stages];

    // Allocate memory for scalars:
    real t;
    real ti;
    real t_inner;
    real hmin = 0.0;
    real trunc_err_max;
    real trunc_err_allowed;
    real delta;
    int step = 0;

    // Initialize:
    t = t0;
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    while (t < tf)
    {
        // printf("\n\nBegin time step t = %f | h = %f\n", t, h);
        hmin = 16.0 * (_dynorb_eps(t)); // Update hmin
        ti = t;
        _dynorb_rvcopy(sys_size, yy, yyi_);
        // Evaluate time derivatives at stage i:
        for (int i = 0; i < n_stages; ++i)
        {
            // Stage i: t_inner and yy_inner
            t_inner = ti + (a[0] * h);
            _dynorb_rvcopy(sys_size, yyi_, yy_inner); // yy_inner = yyi_
            for (int j = 0; j < i; ++j)
            {
                _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, n_stages, i, j)), &FF_[sys_size * j], yy_inner); // yy_inner = yyi_ + b(2,1)*ff1_
            }
            //
            odeFunction(t_inner, yy_inner, odeParams, &FF_[sys_size * i]);
        }

        // Truncation Error Vector:
        _dynorb_rvcopy(sys_size, &FF_[0], ee);
        _dynorb_rscal(sys_size, (c[0] - c_star[0]), ee);
        for (int i = 1; i < n_stages; ++i)
        {
            _dynorb_raxpy(sys_size, (c[i] - c_star[i]), &FF_[sys_size * i], ee);
        }
        _dynorb_rscal(sys_size, h, ee);
        // Maximum Component of Truncation Error:
        trunc_err_max = _dynorb_max_abs_component(sys_size, ee);

        // Allowable Truncation Error:
        trunc_err_allowed = _dynorb_max_abs_component(sys_size, yy);
        trunc_err_allowed = tol * _dynorb_rmax(trunc_err_allowed, 1.0);

        // If truncation error is out bound --> update h:
        if (trunc_err_max > trunc_err_allowed)
        {
            printf("Maximum truncation error = %f\n", trunc_err_max);
            printf("Allowed Truncation error = %f\n", trunc_err_allowed);
            // Fractional Change in step-size:
            delta = (real)pow(trunc_err_allowed / (trunc_err_max + DBL_EPSILON), (1.0 / 5.0)); // COULD GIVE PRROBLEMS float/double
            delta *= 0.8;
            // printf("Fractional change in step-size: delta = %f\n", delta);
            // Update time step :
            // printf("Step size before: h = %f\n", h);
            h = _dynorb_rmin(delta * h, 4.0 * h);
            // printf("Step size after: h = %f\n", h);
            if (h < hmin)
            {
                printf("\n\n Warning: Step size (h = %.10f) fell below\n", h);
                printf("its minimum allowable value (hmin = %.10f) at time %f.\n\n", hmin, t);
                // break;
            }
        }
        else // Update solution
        {
            // Safety check:
            h = _dynorb_rmin(h, fabs(tf - t));
            // Update Time
            t += h;

            // Compute updated solution
            _dynorb_rvcopy(sys_size, &FF_[0], yy);
            _dynorb_rscal(sys_size, c[0], yy);
            for (int i = 1; i < n_stages; ++i)
            {
                _dynorb_raxpy(sys_size, c[i], &FF_[sys_size * i], yy);
                // _dynorb_rmprint(FF_, sys_size*n_stages, 1);
            }
            _dynorb_rscal(sys_size, h, yy);
            _dynorb_raxpy(sys_size, 1.0, yyi_, yy);

            // Realloc if needed:
            if (step >= n_steps)
            {
                // Double the current number of steps (increase size)
                n_steps += 1; // n_steps *= 2;
                // Update Structure:
                solver_configuration->n_steps = n_steps;
                // Reallocate memory for solution arrays
                tt = (real *)realloc(tt, n_steps * sizeof(real));
                YY_t = (real *)realloc(YY_t, n_steps * sys_size * sizeof(real));
            }
            // Store the results

            tt[step] = t; // Store the current time
            // printf("\ntt[step] = %.3f | t = %.3f | ptr_tt = %p\n", tt[step], t, (void *)tt);
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]); // Store the current state
            // Increase step:
            ++step;
            // printf("End time step t = %f | h = %f\n", t, h);
        }
    }
}

/* ====================================================================================================
 * DYNORB ORBITAL MECHANICS FUNCTIONS:
 * ==================================================================================================== */

void _dynorb_twoBodyAbsFun(const real t, const real *yy, const void *params, real *ff)
{
    // Parameters:
    (void)t; // Useless, just to compile
    _dynorb_twoBodyAbsParams *Params = (_dynorb_twoBodyAbsParams *)params;
    real m1 = Params->m1;
    real m2 = Params->m2;
    real G = Params->G;
    // Position:
    real RR1_[3] = {yy[0], yy[1], yy[2]};
    real RR2_[3] = {yy[3], yy[4], yy[5]};
    // Velocity:
    real VV1_[3] = {yy[6], yy[7], yy[8]};
    real VV2_[3] = {yy[9], yy[10], yy[11]};
    // Acceleration:
    real AA1_[3];
    real AA2_[3];
    // Compute Acceleration:
    real RR_diff[3];
    _dynorb_rvcopy(3, RR2_, RR_diff);
    _dynorb_raxpy(3, -1.0, RR1_, RR_diff);
    real r = _dynorb_rnrm2(3, RR_diff);
    _dynorb_rscal(3, (1 / (r * r * r)), RR_diff); // (R2-R1)/r^3
    _dynorb_rvcopy(3, RR_diff, AA1_);             // Compute AA1_ and AA2_
    _dynorb_rvcopy(3, RR_diff, AA2_);             // Compute AA1_ and AA2_
    _dynorb_rscal(3, (G * m2), AA1_);             // Compute AA1_ and AA2_
    _dynorb_rscal(3, (-1.0 * G * m1), AA2_);      // Compute AA1_ and AA2_
    // Update State Derivatives:
    _dynorb_rvcopy(3, VV1_, &ff[0]);
    _dynorb_rvcopy(3, VV2_, &ff[3]);
    _dynorb_rvcopy(3, AA1_, &ff[6]);
    _dynorb_rvcopy(3, AA2_, &ff[9]);
}

void _dynorb_twoBodyRelFun(const real t, const real *yy, const void *params, real *ff)
{
    // Parameters:
    (void)t; // Useless, just to compile
    _dynorb_twoBodyRelParams *Params = (_dynorb_twoBodyRelParams *)params;
    // real m = Params->m;
    real mu = Params->mu;
    // Position:
    real rr_[3] = {yy[0], yy[1], yy[2]};
    // Velocity:
    real vv_[3] = {yy[3], yy[4], yy[5]};
    // Acceleration:
    real aa_[3];
    // Compute Acceleration:
    _dynorb_rvcopy(3, rr_, aa_);
    real r = _dynorb_rnrm2(3, rr_);
    _dynorb_rscal(3, -1.0 * (mu / (r * r * r)), aa_); // (R2-R1)/r^3
    // Update State Derivatives:
    _dynorb_rvcopy(3, vv_, &ff[0]);
    _dynorb_rvcopy(3, aa_, &ff[3]);
}

void _dynorb_threeBodyRestrictFun(const real t, const real *yy, const void *params, real *ff)
{ // Computes the components of the relative acceleration for the restricted 3 - body problem,
    // using Equations 2.192a and 2.192b (pag. 126)

    // Parameters:
    (void)t;
    _dynorb_threeBodyRestrictParams *p = (_dynorb_threeBodyRestrictParams *)params;
    // State Components:
    real x = yy[0];
    real y = yy[1];
    real vx = yy[2];
    real vy = yy[3];
    // Distance from Body 1 and Body 2:
    real r1 = sqrt((pow((x + (p->pi_2 * p->r12)), 2.0)) + (pow(y, 2.0)));
    real r2 = sqrt((pow((x - p->pi_1 * p->r12), 2.0)) + (pow(y, 2.0)));
    // Common factors:
    real _2W = 2.0 * p->W;
    real _W2 = pow(p->W, 2.0);
    real _r1_3 = pow(r1, 3.0);
    real _r2_3 = pow(r2, 3.0);

    // State derivatives of restricted 3-body problem: ff = dyy/dt
    ff[0] = yy[2];
    ff[1] = yy[3];
    ff[2] = _2W * vy + _W2 * x - (p->mu1 * (x - p->x1) / _r1_3) - (p->mu2 * (x - p->x2) / _r2_3);
    ff[3] = -1.0 * _2W * vx + _W2 * y - ((p->mu1 / _r1_3) + (p->mu2 / _r2_3)) * y;
}

real _dynorb_rkeplerH(real e, real M_hyperbola, real tol)
{
    // Initial estimate (rough):
    real F = M_hyperbola;
    real ratio = 1.0;
    while (fabs(ratio) > tol)
    {
        // Update ratio:
        ratio = (e * sinh(F) - F - M_hyperbola) / (e * cosh(F) - 1);
        // Update F:
        F -= ratio;
    }
    return F;
}

real _dynorb_rkeplerE(real e, real M_ellipse, real tol)
{
    // Initial estimate of the root (Prussing & Conway):
    real E;
    if (M_ellipse < _dynorb_PI)
    {
        E = M_ellipse + 0.5 * e;
    }
    else
    {
        E = M_ellipse - 0.5 * e;
    }
    real ratio = 1.0;
    while (fabs(ratio) > tol)
    {
        // Update ratio:
        ratio = (E - e * sin(E) - M_ellipse) / (1 - e * cos(E));
        // Update E:
        E -= ratio;
    }
    return E;
}

void _dynorb_rlagrangeCoefficientsFrom_yy0_Dth(const real mu, const real Dth, const real *yy0, real *fgdfdg)
{
    // Norms
    real r0 = _dynorb_rnrm2(3, &yy0[0]);
    real v0 = _dynorb_rnrm2(3, &yy0[3]);
    // Radial Component of velocity:
    real vr0 = _dynorb_rdot(3, &yy0[0], &yy0[3]) / r0;
    // Magnitude of (constant) angular momentum:
    real h0 = r0 * sqrt((v0 * v0) - (vr0 * vr0));
    real s = sin(Dth);
    real c = cos(Dth);
    // Equation 2.152 :
    real r = ((h0 * h0) / mu) / ((1 + ((((h0 * h0) / (mu * r0)) - 1) * c) - ((h0 * vr0 * s) / mu)));
    // Equations 2.158a & b :
    fgdfdg[0] = 1 - (mu * r * (1 - c) / (h0 * h0));
    fgdfdg[1] = r * r0 * s / h0;
    // Equations 2.158c & d :
    fgdfdg[2] = (mu / h0) * (vr0 / h0 * (1 - c) - s / r0);
    fgdfdg[3] = 1 - mu * r0 / (h0 * h0) * (1 - c);
}

void _dynorb_yy_From_yy0_Dth(const real mu, const real Dth, const real *yy0, real *yy)
{
    real fgdfdg[4];
    _dynorb_rlagrangeCoefficientsFrom_yy0_Dth(mu, Dth, yy0, fgdfdg);
    // Computation of rr:
    _dynorb_rvcopy(3, &yy0[0], &yy[0]);
    _dynorb_rscal(3, fgdfdg[0], &yy[0]);
    _dynorb_raxpy(3, fgdfdg[1], &yy0[3], &yy[0]);
    // Computation of vv:
    _dynorb_rvcopy(3, &yy0[0], &yy[3]);
    _dynorb_rscal(3, fgdfdg[2], &yy[3]);
    _dynorb_raxpy(3, fgdfdg[3], &yy0[3], &yy[3]);
}

real _dynorb_rkeplerUniversal(const real mu, const real Dt, const real r0, const real vr0, const real alpha, const real tol, const int max_iter)
{
    // Startig value for Universal Variable:
    real sqrt_mu = sqrt(mu);
    real chi = sqrt_mu * fabs(alpha) * Dt;
    // Initialize ratio and number of iterations:
    real ratio = 1.0;
    int n_iter = 0;
    // printf("\n_dynorb_rkeplerUniversal: Iteration [%d] | ratio = %g | Chi = %g\n", n_iter, ratio, chi);

    while (fabs(ratio) > tol && n_iter < max_iter)
    {
        // Temporary variables:
        real A = (r0 * vr0) / sqrt_mu;
        real B = 1.0 - (alpha * r0);
        real chi2 = chi * chi;
        real chi3 = pow(chi, 3.0);
        real z = alpha * chi2;
        real S = _dynorb_rstumpffS(z);
        real C = _dynorb_rstumpffC(z);
        // Update ratio
        ratio = ((A * chi2 * C) + (B * chi3 * S) + (r0 * chi) - (sqrt_mu * Dt)) /
                ((A * chi * (1 - alpha * chi2 * S)) + (B * chi2 * C) + r0);
        // Update Universal Variable:
        chi -= ratio;
        // Update number of iter:
        ++n_iter;
        // printf("_dynorb_rkeplerUniversal: Iteration [%d] | ratio = %g | Chi = %g\n", n_iter, ratio, chi);
    }
    if (n_iter > max_iter)
    {
        printf("\n");
        printf("_dynorb_rkeplerUniversal: Number of iterations (n_iter = %d) reached limit: %d\n", n_iter, max_iter);
    }
    return chi;
}

void _dynorb_rfgFromChi(const real mu, const real Dt, const real chi, const real r0, const real alpha, real *fg)
{ // Computes Lagrange coefficients from Universal Variable
    // Temporary variables:
    real chi2 = chi * chi;
    real chi3 = pow(chi, 3.0);
    real sqrt_mu = sqrt(mu);
    // Variable z:
    real z = alpha * chi2;
    // Stumpff Functions:
    real S = _dynorb_rstumpffS(z);
    real C = _dynorb_rstumpffC(z);
    // Lagarange coefficients (Formula 3.69, pag.181):
    fg[0] = 1.0 - ((chi2 / r0) * C);
    fg[1] = Dt - ((chi3 * S) / sqrt_mu);
}
void _dynorb_rdfdgFromChi(const real mu, const real chi, const real r, const real r0, const real alpha, real *dfdg)
{ // Computes time-derivatives of Lagrange coefficients from Universal Variable
    // Temporary variables:
    real chi2 = chi * chi;
    real chi3 = pow(chi, 3.0);
    real sqrt_mu = sqrt(mu);
    // Variable z:
    real z = alpha * chi2;
    // Stumpff Functions:
    real S = _dynorb_rstumpffS(z);
    real C = _dynorb_rstumpffC(z);
    // Lagarange coefficients (Formula 3.69, pag.181):
    dfdg[0] = (sqrt_mu / (r * r0)) * ((alpha * chi3 * S) - chi);
    dfdg[1] = 1.0 - ((chi2 / r) * C);
}

void _dynorb_ryy_From_yy0_Dt(const real mu, const real Dt, const real *yy0, real *yy, const real tol, const int max_iter)
{
    // Vector Norms:
    real r0 = _dynorb_rnrm2(3, &yy0[0]);
    real v0 = _dynorb_rnrm2(3, &yy0[3]);
    // Radial Component of Velocity :
    real vr0 = (_dynorb_rdot(3, &yy0[0], &yy0[3])) / r0;
    // Reciprocal of semi-major axis:
    real alpha = (2.0 / r0) - ((v0 * v0) / mu);
    // Universal Anomaly Chi:
    real chi = _dynorb_rkeplerUniversal(mu, Dt, r0, vr0, alpha, tol, max_iter);
    // Compute f and g:
    real fg[2];
    _dynorb_rfgFromChi(mu, Dt, chi, r0, alpha, fg);
    // Final Position Vector rr = f*rr0 + g*vv0:
    _dynorb_rvcopy(3, &yy0[0], &yy[0]);
    _dynorb_rscal(3, fg[0], &yy[0]);
    _dynorb_raxpy(3, fg[1], &yy0[3], &yy[0]);
    // Norm of rr:
    real r = _dynorb_rnrm2(3, &yy[0]);
    // Compute time-derivatives f and g:
    real dfdg[2];
    _dynorb_rdfdgFromChi(mu, chi, r, r0, alpha, dfdg);
    // Final Velocity Vector vv = df*rr0 + dg*vv0:
    _dynorb_rvcopy(3, &yy0[0], &yy[3]);         // vv<--rr0
    _dynorb_rscal(3, dfdg[0], &yy[3]);          // vv<--df*vv
    _dynorb_raxpy(3, dfdg[1], &yy0[3], &yy[3]); //  vv<--vv + dg*vv0
}
#endif // DYNORB_IMPLEMENTATION
