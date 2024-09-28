/* ****************************************************
 * ****************** DYNORB HEADER *******************
 * **************************************************** */

/* TODO: REMOVE THESE 2 lines: */
#define USE_DOUBLE
#define USE_CBLAS

/* ASSUMPTIONS:
 * 1. Code uses COLUMN MAJOR convention (see CBLAS).
 * 2. Strides other than 1 are not supported in this implementation.
 * 3. Type 'real' can only be float (single precision) or double (double precision).
 */

/* TODO: Short/Medium Term
 * 1. Implement Curtis Algorithms using wrappers of CBLAS and LAPACKE
 * 2. Reduce (completly eliminate when possible) heap usage
 * 3. Respect as much as possible NASA POWER OF TEN
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

/* =====================================
 * SATNDARD C-LIBS:
 * ===================================== */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* =====================================
 * CBLAS SWITCH:
 * ===================================== */
#ifdef USE_CBLAS
#include <cblas.h>
#endif

/* =====================================
 * DOUBLE-FLOAT SWITCH:
 * ===================================== */
#ifdef USE_DOUBLE
typedef double real;
const bool real_is_double = 1;
#elif defined(USE_FLOAT)
typedef float real;
const bool real_is_double = 0;
#endif

/* =====================================
 * DYNORB MACROS:
 * ===================================== */
#define _dynorb_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#define _dynorb_MU_E 398600.44188

/* =====================================
 * DYNORB GLOBAL VARIABLES:
 * ===================================== */

/* =====================================
 * DYNORB CORE FYNCTION TYPES:
 * ===================================== */
typedef real (*_dynorb_nonLinScalFun)(real a, void *funParams);                                 // Non-linear function
typedef void(_dynorb_odeFun)(const real t, const real *yy, const void *odeParams, real *dyydt); // ODE-function

/* =====================================
 * DYNORB CORE STRUCTURES:
 * ===================================== */
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

/* =====================================
 * DYNORB USER SUPPORT STRUCTURES:
 * ===================================== */
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

/* =====================================
 * DYNORB UTILS FUNCTIONS:
 * ===================================== */

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

/* ===========================================
 * DYNORB BASIC LINEAR ALGEBRA SUBROUTINES:
 * =========================================== */
// Thin wrappers of CBLAS if USE_CBLAS. Else custum implementation.
real _dynorb_rdot(const int n, const real *xx, const real *yy);
real _dynorb_rnrm2(const int n, const real *xx);
void _dynorb_rscal(const int n, const real alpha, real *xx);
void _dynorb_raxpy(const int n, const real alpha, const real *xx, real *yy);
void _dynorb_rgemv(const bool TransposeA, const int m, const int n, const real alpha, const real *A, const real *xx, const real beta, real *yy);
void _dynorb_rvcopy(const int n, const real *src_xx, real *dst_yy);
void _dynorb_rmcopy(const real *src_A, const int src_m, const int src_n, const int src_start_i, const int src_start_j, const int cpy_m, const int cpy_n, real *dst_B, const int dst_m, const int dst_n, const int dst_start_i, const int dst_start_j);

/* =====================================
 * DYNORB USER SUPPORT FUNCTIONS:
 * ===================================== */

/**
 * @brief Computes the state derivatives for the restricted 3-body problem in a rotating frame.
 *
 * This function calculates the relative acceleration and state derivatives for a test particle
 * in the restricted 3-body problem, where two massive bodies (mu1, mu2) influence the motion of
 * a massless third body. The equations are derived from Equations 2.192a and 2.192b (Curtis, pg. 126).
 *
 * @param t         Time variable (unused in calculations).
 * @param yy        Pointer to an array of state variables (4 elements):
 *                  - Position in rotating frame: yy[0] = x, yy[1] = y
 *                  - Velocity in rotating frame: yy[2] = vx, yy[3] = vy
 * @param params    Pointer to a structure containing the problem parameters:
 *                  - W: Angular velocity of the rotating frame
 *                  - r12: Distance between the two main bodies
 *                  - pi_1, pi_2: Proportional distances of the main bodies
 *                  - mu1, mu2: Gravitational parameters of the main bodies
 *                  - x1, x2: X-coordinates of the main bodies in the rotating frame
 * @param ff        Pointer to an array where the computed derivatives will be stored (4 elements):
 *                  - ff[0] = vx (derivative of position x)
 *                  - ff[1] = vy (derivative of position y)
 *                  - ff[2] = ax (derivative of velocity vx)
 *                  - ff[3] = ay (derivative of velocity vy)
 *
 * This function computes the derivatives of both the position and velocity of the test particle
 * under the influence of the gravitational forces from two primary bodies in a rotating reference frame.
 */
void _dynorb_threeBodyRestrictFun(const real t, const real *yy, const void *params, real *ff);

/**
 * @brief Propagates the initial state vector \p yy0 by a true anomaly change \p Dth
 *        to compute the final state vector \p yy using Lagrange's equations.
 *
 * This function computes the new position and velocity vectors by using the
 * Lagrange coefficients, which are obtained from the initial state and the
 * true anomaly difference.
 *
 * @param[in]  mu   Gravitational parameter (standard gravitational parameter) of the central body.
 * @param[in]  Dth  Change in true anomaly in radians.
 * @param[in]  yy0  Initial state vector (position and velocity), an array of size 6:
 *                    - First 3 elements are the position vector.
 *                    - Last 3 elements are the velocity vector.
 * @param[out] yy   Final state vector (position and velocity), an array of size 6:
 *                    - First 3 elements are the updated position vector.
 *                    - Last 3 elements are the updated velocity vector.
 *
 * This function internally uses Lagrange's coefficients to compute the final state vector.
 */
void _dynorb_yy_From_yy0_Dth(const real mu, const real Dth, const real *yy0, real *yy);

/**
 * @brief Finds the root of a nonlinear scalar function using the Bisection Method.
 *        // TODO: WORKS ONLY WITH DOUBLES for now
 *
 * This function implements the Bisection Method to approximate the root of a scalar function
 * within a given interval [a, b]. It iteratively halves the interval based on the sign of the
 * function evaluated at the midpoint until the desired tolerance is achieved.
 *
 * @param function  Pointer to the function to be solved. The function should take two arguments:
 *                  - A real-valued input (the variable of the function).
 *                  - A pointer to additional parameters (passed through funParams).
 * @param funParams Pointer to any additional parameters that need to be passed to the function.
 * @param a         The lower bound of the interval.
 * @param b         The upper bound of the interval.
 * @param tol       The tolerance for the method. The algorithm will stop when the interval size
 *                  is smaller than this tolerance.
 * @return          The midpoint of the final interval, which is the approximation of the root.
 *
 * This function uses the Bisection Method to find a root of the nonlinear function within
 * the specified interval [a, b]. The method works by narrowing down the interval in which the
 * root is located based on the signs of the function at the boundaries. The process continues
 * until the interval is smaller than the tolerance.
 *
 * The function assumes that the function values at the initial points a and b have opposite signs
 * (i.e., f(a) * f(b) < 0), ensuring that a root exists in the interval.
 */
real _dynorb_bisect(_dynorb_nonLinScalFun function, void *funParams, real a, real b, real tol);

/**
 * @brief Computes the Lagrange coefficients for orbital propagation given the initial state and true anomaly change.
 *
 * This function computes the Lagrange coefficients \( f \), \( g \), \( f_dot \), and \( g_dot \)
 * used to propagate the orbital state using Lagrange's equations.
 *
 * @param[in]  mu      Gravitational parameter (standard gravitational parameter) of the central body.
 * @param[in]  Dth     Change in true anomaly in radians.
 * @param[in]  yy0     Initial state vector (position and velocity), an array of size 6:
 *                      - First 3 elements are the position vector.
 *                      - Last 3 elements are the velocity vector.
 * @param[out] fgdfdg  Array of size 4 to store the Lagrange coefficients:
 *                      - fgdfdg[0] = \( f \)
 *                      - fgdfdg[1] = \( g \)
 *                      - fgdfdg[2] = \( \dot{f} \)
 *                      - fgdfdg[3] = \( \dot{g} \)
 *
 * This function calculates the radial and tangential components of velocity, as well as the constant angular momentum,
 * and applies the appropriate equations to compute the Lagrange coefficients.
 */
void _dynorb_LagrangeFunctionsFrom_yy0_Dth(const real mu, const real Dth, const real *yy0, real *fgdfdg);

/**
 * @brief Computes the RELATIVE derivatives (in NON-ROTATING frame) for a secondary body orbiting a main attractor.
 *
 * This function calculates the RELATIVE acceleration and state derivatives the secondary body
 * under the gravitational attraction of the main one. The quantities are relative (to main body) and
 * expressed in a NON-ROTATING frame centered in the main body (CoM). For Earth orbits, this is the ECI.
 * The function takes the current state of the system, parameters, and outputs the derivatives of the
 * state variables. See Curtis,  Chapter 2, pag.69: Algorithm 2.1
 *
 * @param t          Time variable (unused in calculations).
 * @param yy        Pointer to an array of state variables (12 elements):
 *                  - Relative Position of secondary body: yy[0], yy[1], yy[2]
 *                  - Relative Velocity of secondary body: yy[3], yy[4], yy[5]
 * @param params    Pointer to a structure containing gravitational parameters:
 *                  - m: Mass of secondary body
 *                  - mu: Gravitational parameter of main attractor
 * @param ff        Pointer to an array where the computed derivatives will be stored:
 *                  - Derivative of secondary body velocity: ff[0], ff[1], ff[2]
 *                  - Derivative of secondary body acceleration: ff[6], ff[7], ff[8]
 *
 * The function computes the acceleration for each body based on their positions
 * and outputs the derivatives of both position and velocity.
 */
void _dynorb_twoBodyRelFun(const real t, const real *yy, const void *params, real *ff);

/**
 * @brief Computes the ABSOLUTE derivatives for a two-body dynamical system.
 *
 * This function calculates the ABSOLUTE acceleration and state derivatives for two bodies
 * under the influence of their mutual gravitational attraction. Quantities are expressed in the
 * Inertial frame coordinate system.
 * The function takes the current state of the system, parameters, and outputs the derivatives of the
 * state variables. See Curtis,  Chapter 2, pag.64: Algorithm 2.1
 *
 * @param t          Time variable (unused in calculations).
 * @param yy        Pointer to an array of state variables (12 elements):
 *                  - Position of body 1: yy[0], yy[1], yy[2]
 *                  - Position of body 2: yy[3], yy[4], yy[5]
 *                  - Velocity of body 1: yy[6], yy[7], yy[8]
 *                  - Velocity of body 2: yy[9], yy[10], yy[11]
 * @param params    Pointer to a structure containing gravitational parameters:
 *                  - G: Gravitational constant
 *                  - m1: Mass of body 1
 *                  - m2: Mass of body 2
 * @param ff        Pointer to an array where the computed derivatives will be stored:
 *                  - Derivative of body 1 velocity: ff[0], ff[1], ff[2]
 *                  - Derivative of body 2 velocity: ff[3], ff[4], ff[5]
 *                  - Derivative of body 1 acceleration: ff[6], ff[7], ff[8]
 *                  - Derivative of body 2 acceleration: ff[9], ff[10], ff[11]
 *
 * The function computes the acceleration for each body based on their positions
 * and outputs the derivatives of both position and velocity.
 */
void _dynorb_twoBodyAbsFun(const real t, const real *yy, const void *params, real *ff);

/* =====================================
 * DYNORB USER SUPPORT FUNCTIONS:
 * ===================================== */

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
void _dynorb_rcol(const real *A, int m, int j, real *column);

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
void _dynorb_rrow(const real *A, int m, int n, int i, real *row);

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
static inline real _dynorb_rel(const real *A, int m, int i, int j);

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
 * @param[in] tf Final time of the simulation.
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
                               const int sys_size, const real t0, const real tf, const real h);

/**
 * @brief Configures an ODE (Ordinary Differential Equation) system assuming stack memory allocation.
 *
 * This function sets up an ODE system for solving using stack memory allocation. It calculates the
 * number of steps required and provides guidelines for manually allocating memory on the stack within
 * the user's code. The function does not perform memory allocation itself but instead assumes that the
 * user will handle it. To manually allocate memory on the stack for the ODE system, include the following lines in your code:
 *
 * - real tt[solver_configuration.n_steps];
 *
 * - real YY_t[solver_configuration.n_steps * ode_system.sys_size];
 *
 * - ode_system.tt = tt;
 *
 * - ode_system.YY_t = YY_t;
 *
 * @param[in] ode_system Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to the solver configuration structure.
 * @param[in] ode_function Pointer to the ODE function that needs to be solved.
 * @param[in] odeParams Pointer to any additional parameters needed by the ODE function.
 * @param[in] yy0 Pointer to the initial conditions of the system.
 * @param[in] sys_size The size of the system (i.e., the number of state variables).
 * @param[in] t0 Initial time of the simulation.
 * @param[in] tf Final time of the simulation.
 * @param[in] h The time step size.
 *
 * This function assumes the use of a column-major ordering for matrix data.
 * Stack memory allocation is expected to be handled by the user, as demonstrated in the code block above.
 *
 * @return[out] void
 */
void _dynorb_configure_static(_dynorb_odeSys *ode_system, _dynorb_solverConf *solver_configuration,
                              _dynorb_odeFun *ode_function, void *odeParams, real *yy0,
                              const int sys_size, const real t0, const real tf, const real h);

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
void _dynorb_rrk1(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

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
void _dynorb_rrk2(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

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
void _dynorb_rrk3(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

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
void _dynorb_rrk4(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration);

/**
 *
 * @brief Performs Heun Predictor-Corrector numerical integration for ODE systems.
 *
 * This function integrates an ODE system using Heun Predictor-Corrector method.
 *
 * @param[inout] sys Pointer to the structure defining the ODE system.
 * @param[in] solver_configuration Pointer to structure defining Solver Configuration.
 * @param[in] tol Tolerance
 * @param[in] max_iter Maximum Number of Prediction-Correction Iterations
 *
 * This function uses a column-major convention.
 *
 * @return Void.
 *
 */
void _dynorb_rheun(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration, const real tol, const int max_iter);

#endif // DYNORB_H_
