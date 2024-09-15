/* ****************************************************
 * ****************** DYNORB HEADER *******************
 * **************************************************** */

/* ASSUMPTIONS:
 * 1. Code uses COLUMN MAJOR convention (see CBLAS).
 * 2. Strides other than 1 are not supported in this implementation.
 * 3. Type 'real' can only be float (single precision) or double (double precision).
 */

/* TODO: Short/Medium Term
 * 1. Implement Curtis Algorithms usinng wrappers of CBLAS and LAPACKE
 * 2. Reduce (completly eliminate when possible) heap usage
 * 3. Fullfill as much as possible NASA POWER OF TEN
 */

/* TODO: Long Term
 * 1. Implement Custom subroutines to substitute CBLAS and LAPACKE
 * 2. Completly eliminate heap usage (if impossible, remove the function from library)
 * 3. Enforce NASA POWER OF TEN
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

/* ODE FUNCTION: */
typedef void(_dynorb_odeFun)(const real t, const real *yy, const void *odeParams, real *dyydt);

/* ODE SYSTEM AND SOLVER: */
typedef struct
{
    _dynorb_odeFun *odeFunction; // Pointer to the ODE function
    void *odeParams;             // Pointer to the parameters (structure) for the ODE function
    int sys_size;                // Size of the system (number of equations)
    real t0;                     // Initial time
    real t1;                     // Final time
    real *yy0;                   // Pointer to the initial state array
    real *tt;                    // Time Steps of the Solution
    real *YY_t;                  // Solution Array
} _dynorb_odeSys;

typedef struct
{
    real t0;     // Initial time
    real t1;     // Final time
    real h;      // Initial Step Size
    int n_steps; // Initial Number of steps
} _dynorb_solverConf;

/* UTILITY FUNCTIONS: */
int _dynorb_min(int a, int b);
int _dynorb_max(int a, int b);
real _dynorb_rmin(real a, real b);
real _dynorb_rmax(real a, real b);
void _dynorb_rmprint(const real *A, const int m, const int n);
void _dynorb_rvfill(real *xx, const int n, const real a);

/* BASIC LINEAR ALGEBRA SUBROUTINES: */
// Thin wrappers of CBLAS if USE_CBLAS. Else custum implementation.
real _dynorb_rdot(const int n, const real *xx, const real *yy);
real _dynorb_rnrm2(const int n, const real *xx);
void _dynorb_rscal(const int n, const real alpha, real *xx);
void _dynorb_raxpy(const int n, const real alpha, const real *xx, real *yy);
void _dynorb_rgemv(const bool TransposeA, const int m, const int n, const real alpha, const real *A, const real *xx, const real beta, real *yy);
void _dynorb_rvcopy(const int n, const real *src_xx, real *dst_yy);
void _dynorb_rmcopy(const real *src_A, const int src_m, const int src_n, const int src_start_i, const int src_start_j, const int cpy_m, const int cpy_n, real *dst_B, const int dst_m, const int dst_n, const int dst_start_i, const int dst_start_j);

/* CORE FUNCTIONS: */

/**
 * @brief Extracts a column from a matrix stored in column-major order.
 *
 * This function extracts the column at index `j` from the matrix `A` (stored in
 * column-major order) and stores it in the array `column`.
 *
 * @param[in] A Pointer to the matrix data (stored in column-major order).
 * @param[in] m Number of rows in the matrix.
 * @param[in] j Column index to extract.
 * @param[out] column Array to store the extracted column (size should be `m`).
 *
 * @return[out] void
 */
void _dynorb_col(const real *A, int m, int j, real *column);

/**
 * @brief Extracts a row from a matrix stored in column-major order.
 *
 * This function extracts the row at index `i` from the matrix `A` (stored in
 * column-major order) and stores it in the array `row`.
 *
 * @param[in] A Pointer to the matrix data (stored in column-major order).
 * @param[in] m Number of rows in the matrix.
 * @param[in] n Number of columns in the matrix.
 * @param[in] i Row index to extract.
 * @param[out] row Array to store the extracted row (size should be `n`).
 *
 * @return[out] void
 */
void _dynorb_row(const real *A, int m, int n, int i, real *row);

/**
 * @brief Accesses an element from a matrix stored in column-major order.
 *
 * This function retrieves the element at position (i, j) from a matrix `M`
 * stored in column-major order.
 *
 * @param[in] M Pointer to the matrix data (assumed to be stored in column-major order).
 * @param[in] m Number of rows in the matrix.
 * @param[in] i Row index of the element.
 * @param[in] j Column index of the element.
 *
 * @return[out] real The value at the (i, j) position in the matrix.
 */
static inline real _dynorb_el(const real *A, int m, int i, int j);

/**
 * @brief Frees dynamically allocated memory in the ODE system.
 *
 * This function releases the memory allocated for the time steps array (`tt`)
 * and the solution array (`YY_t`) in the ODE system.
 *
 * @param[in] ode_system Pointer to the structure defining the ODE system.
 *
 * @return[out] void
 */
void _dynorb_free(_dynorb_odeSys *ode_system);

/**
 * @brief Configures an ODE system using dynamic memory allocation.
 *
 * This function initializes the ODE system for solving with dynamic memory
 * allocation for time steps (`tt`) and solution array (`YY_t`). It computes
 * the number of steps, allocates memory, and fills the ODE system and solver
 * configuration structures. It also checks for memory allocation failures.
 *
 * @param[in] ode_system Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to the solver configuration structure.
 * @param[in] ode_function Pointer to the ODE function to be solved.
 * @param[in] odeParams Pointer to any additional parameters for the ODE function.
 * @param[in] yy0 Pointer to the initial conditions of the system.
 * @param[in] sys_size Size of the system (number of state variables).
 * @param[in] t0 Initial time of the simulation.
 * @param[in] t1 Final time of the simulation.
 * @param[in] h Time step size.
 *
 * This function uses a column-major convention.
 * This function assumes that the user will manually free the allocated memory using
 * `_dynorb_free()`.
 *
 * @return[out] void
 */
void _dynorb_configure_dynamic(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration,
                               _dynorb_odeFun *ode_function, void *odeParams, real *yy0,
                               const int sys_size, const real t0, const real t1, const real h);

/**
 * @brief Configures an ODE system assuming stack memory allocation.
 *
 * This function initializes the ODE system for solving with stack memory
 * allocation. It computes the number of steps and provides instructions
 * for manually allocating stack memory in the user's code.
 *
 * @param[in] ode_system Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to the solver configuration structure.
 * @param[in] ode_function Pointer to the ODE function to be solved.
 * @param[in] odeParams Pointer to any additional parameters for the ODE function.
 * @param[in] yy0 Pointer to the initial conditions of the system.
 * @param[in] sys_size Size of the system (number of state variables).
 * @param[in] t0 Initial time of the simulation.
 * @param[in] t1 Final time of the simulation.
 * @param[in] h Time step size.
 *
 * This function uses a column-major convention.
 * This function assumes that the user will allocate stack memory in their code.
 *
 * @return[out] void
 */
void _dynorb_configure_static(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration,
                              _dynorb_odeFun *ode_function, void *odeParams, real *yy0,
                              const int sys_size, const real t0, const real t1, const real h);

/**
 *
 * @brief Performs Euler (RK1) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Euler (Runge-Kutta order 1) method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to structure defining Solver Configuration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk1(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

/**
 *
 * @brief Performs Heun (RK2) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Heun (Runge-Kutta order 2) method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to structure defining Solver Configuration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk2(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

/**
 *
 * @brief Performs Heun (RK3) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Runge-Kutta order 3 method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to structure defining Solver Configuration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk3(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

/**
 *
 * @brief Performs Heun (RK4) numerical integration for ODE systems.
 *
 * This function integrates an ODE system using a specified order of the Runge-Kutta order 4 method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to structure defining Solver Configuration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rk4(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

/**
 *
 * @brief Performs Heun Predictor-Corrector numerical integration for ODE systems.
 *
 * This function integrates an ODE system using Heun Predictor-Corrector method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to structure defining Solver Configuration.
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_heun_(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

#endif // DYNORB_H_

/*
 * **************************************************** *
 * ************* DYNORNB IMPLEMENTATION ************* *
 * **************************************************** *
 */

#ifdef DYNORB_IMPLEMENTATION

/* BASIC LINEAR ALGEBRA SUBROUTINES: */
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

#endif

#ifndef USE_CBLAS // --> BLAS-like Custom Implementation
/* TODO: custom implementation-> */
void IMPLEMENT_ALTERNATIVE_TO_CBLAS_FUNCTIONS(void)
{
    printf("I HAVE TO IMPLEMENT ALTERNATIVE TO CBLAS FUNCTIONS.\n");
}
#endif

/* UTILITY FUNCTIONS: */
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
            printf("%.4f", _dynorb_el(A, m, i, j));
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

/* CORE FUNCTIONS: */

void _dynorb_col(const real *A, int m, int j, real *column)
{
    for (int i = 0; i < m; i++)
    {
        column[i] = _dynorb_el(A, m, i, j); // Accessing element A[i, j]
    }
}

void _dynorb_row(const real *A, int m, int n, int i, real *row)
{
    for (int j = 0; j < n; j++)
    {
        row[j] = _dynorb_el(A, m, i, j); // Accessing element A[i, j]
    }
}

static inline real _dynorb_el(const real *A, int m, int i, int j)
{
    return A[(j * m) + i];
}

void _dynorb_free(_dynorb_odeSys *ode_system)
{
    free(ode_system->tt);
    free(ode_system->YY_t);
}

void _dynorb_configure_dynamic(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration,
                               _dynorb_odeFun *ode_function, void *odeParams, real *yy0,
                               const int sys_size, const real t0, const real t1, const real h)
{
    // Compute (initial) number of steps:
    int n_steps = (int)((((t1 - t0)) / h) + 1.0);

    // Allocate memory for the time steps and solution array
    ode_system->tt = (real *)malloc(n_steps * sizeof(real));
    ode_system->YY_t = (real *)malloc(n_steps * sys_size * sizeof(real));

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
    ode_system->t1 = t1;
    ode_system->yy0 = yy0;

    // Copy values in solver_configuration structure:
    solver_configuration->t0 = t0;
    solver_configuration->t1 = t1;
    solver_configuration->h = h;
    solver_configuration->n_steps = n_steps;

    // Print:
    printf("\nTime Step Size: <h = %f [s]>\n", solver_configuration->h);
    printf("Number of Steps: <solv_conf.n_steps = %d>\n", solver_configuration->n_steps);
}

void _dynorb_configure_static(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration,
                              _dynorb_odeFun *ode_function, void *odeParams, real *yy0,
                              const int sys_size, const real t0, const real t1, const real h)
{
    // Compute (initial) number of steps:
    int n_steps = (int)((((t1 - t0)) / h) + 1.0);

    // Print stack memory allocation instructions
    ode_system->tt = NULL;
    ode_system->YY_t = NULL;
    printf("\nFor Stack Memory Allocation add these lines to your code:\n");
    printf("    real tt[solver_configuration.n_steps];\n");
    printf("    real YY_t[solver_configuration.n_steps * ode_system.sys_size];\n");
    printf("    ode_system.tt = tt;\n");
    printf("    ode_system.YY_t = YY_t;\n\n");

    // Copy values in ode_system structure:
    ode_system->odeFunction = ode_function;
    ode_system->odeParams = odeParams;
    ode_system->sys_size = sys_size;
    ode_system->t0 = t0;
    ode_system->t1 = t1;
    ode_system->yy0 = yy0;

    // Copy values in solver_configuration structure:
    solver_configuration->t0 = t0;
    solver_configuration->t1 = t1;
    solver_configuration->h = h;
    solver_configuration->n_steps = n_steps;

    // Print:
    printf("\nTime Step Size: <h = %f [s]>\n", solver_configuration->h);
    printf("Number of Steps: <solv_conf.n_steps = %d>\n", solver_configuration->n_steps);
}

void _dynorb_rk1(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Open up _dynorb_solvConf [TODO: remove this]
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    //  RK1 (Euler)
    const int n_stages = 1;
    real a = 0.0; // Note: a is [1 x 1]
    // real b = 1.0; // Note: b is [1 x 1]
    real c = 1.0; // Note: c is [1 x 1]

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size];
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real ff[sys_size * n_stages]; // (ff = dyy/dt)

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives in [t, t+h]
        real t_inner = t + a * h;
        odeFunction(t_inner, yy, odeParams, ff);

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
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
            printf("_dynorb_rk1 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        else
        {
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
        }
    }
}

void _dynorb_rk2(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Open up _dynorb_solvConf [TODO: remove this]
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    //  RK2 (Heun)
    const int n_stages = 2;
    real a[2] = {0.0, 0.0}; // Note: a is [1 x 2]
    real b[2] = {0.0, 0.0}; // Note: b is [2 x 1]
    real c[2] = {0.0, 0.0}; // Note: c is [2 x 1]
    a[1] = 1.0;
    b[0] = 0.0;
    b[1] = 1.0;
    c[0] = 0.5;
    c[1] = 0.5;

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size];
    real yy_inner[sys_size];
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real ff[sys_size * n_stages];

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * h;
            _dynorb_rvcopy(sys_size, yy, yy_inner);

            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, odeParams, &ff[i * sys_size]);
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
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
            printf("_dynorb_rk2 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        else
        {
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
        }
    }
}

void _dynorb_rk3(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Open up _dynorb_solvConf [TODO: remove this]
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    //  RK3
    const int n_stages = 3;
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

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size];
    real yy_inner[sys_size];
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real ff[sys_size * n_stages];

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * h;
            _dynorb_rvcopy(sys_size, yy, yy_inner);
            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, odeParams, &ff[i * sys_size]);
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
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
            printf("_dynorb_rk3 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        else
        {
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
        }
    }
}

void _dynorb_rk4(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up _dynorb_odeSys [TODO: remove this]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Open up _dynorb_solvConf [TODO: remove this]
    const real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    //  RK4 (Runge-Kutta)
    const int n_stages = 4;
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
    real yy[sys_size];
    real yy_inner[sys_size];
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Allocate Memory for derivatives at each stage:
    real ff[sys_size * n_stages];

    // Integration Time Instant:
    real t = t0;

    // Numerical Integration:
    for (int step = 0; step < n_steps; ++step)
    {
        // Evaluate Time Derivatives at 'n_stages' points in [t, t+h]
        for (int i = 0; i < n_stages; ++i)
        {
            real t_inner = t + a[i] * h;
            _dynorb_rvcopy(sys_size, yy, yy_inner);
            for (int j = 0; j < i; ++j)
            {
                for (int k = 0; k < sys_size; ++k)
                {
                    yy_inner[k] += h * b[i * (n_stages - 1) + j] * ff[j * sys_size + k];
                }
            }
            odeFunction(t_inner, yy_inner, odeParams, &ff[i * sys_size]);
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
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
            printf("_dynorb_rk4 : Breaking From Loop<t = %f [s]>", t);
            break;
        }
        else
        {
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
            tt[step] = t;
        }
    }
}

void _dynorb_heun_(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration)
{
    // Open up _dynorb_odeSys
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array

    // Open up _dynorb_solvConf [TODO: remove this]
    real h = solver_configuration->h;
    const int n_steps = solver_configuration->n_steps;

    // Tolerance and _dynorb_Max number of steps:
    const real tol = 1.e-6;
    const int rmax_iter = 100;

    // Initialize Algorithm:
    real t = t0;
    real yy[sys_size];
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    // Copy First State and Instant in Output:
    tt[0] = t0;
    _dynorb_rvcopy(sys_size, yy0, YY_t); // Current state at t0 = initial state

    // Allocate Memory for state and derivatives at interval boundaries:
    real yy1_[sys_size];
    real yy2_[sys_size];
    real ff1_[sys_size];
    real ff2_[sys_size];
    real yy2pred_[sys_size];
    real ffavg_[sys_size];

    // Main Loop
    int step = 1;
    while (t < t1 && step < n_steps)
    {
        // Step Size:
        h = _dynorb_rmin(h, t1 - t);

        // Left Boundary of Interval:
        real t1_ = t;
        _dynorb_rvcopy(sys_size, yy, yy1_);
        odeFunction(t1_, yy1_, odeParams, ff1_);

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
            _dynorb_rvcopy(sys_size, yy2_, yy2pred_);
            odeFunction(t2_, yy2pred_, odeParams, ff2_);

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

            // Error calculation // TODO: remove fabs
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
        _dynorb_rvcopy(sys_size, yy2_, yy);
        tt[step] = t;
        _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]);
        step++;
    }
}

#endif // DYNORB_IMPLEMENTATION
