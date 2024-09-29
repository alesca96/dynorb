#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

/*  Example 3.7 (Chapter 3, pag.181):This program uses Algorithm 3.4 and the data
    of Example 3.7 to solve Kepler’s equation.
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

int main(void)
{
    // DATA:
    real mu_E = (real)_dynorb_MU_E;   // Earth grav. parameter [km^3 /s^2]
    real rr0[3] = {7000.0, -12124.0}; // Initial Position vector, ECI frame [km]
    real vv0[3] = {2.6679, 4.6210};   // Initial Velocity vector, ECI frame [km/s]
    real Dt = 1.0 * 60 * 60;          // Time Inteval from (t0, rr0, vv0) [s]

    // Tolerance and Max iterations:
    real tol = 1.e-8;
    int max_iter = 100;

    // EXERCISE 3.7:
    // Initial State:
    real yy0[6];
    _dynorb_rvcopy(3, rr0, &yy0[0]);
    _dynorb_rvcopy(3, vv0, &yy0[3]);
    real yy[6];
    _dynorb_ryy_From_yy0_Dt(mu_E, Dt, yy0, yy, tol, max_iter);

    // Print Info:
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.7: \n");
    printf("-----------------------\n\n");
    printf("Initial Position Vector [km]:\nrr0 = \n");
    _dynorb_rmprint(&yy0[0], 3, 1);
    printf("\n");
    printf("Initial Velocity Vector [km/s]:\nvv0 = \n");
    _dynorb_rmprint(&yy0[3], 3, 1);
    printf("\n");
    printf("Elapsed Time Dt = %g [s]\n", Dt);
    printf("\n");
    printf("Final Position Vector [km]:\nrr = \n");
    _dynorb_rmprint(&yy[0], 3, 1);
    printf("\n");
    printf("Final Velocity Vector [km/s]:\nvv = \n");
    _dynorb_rmprint(&yy[3], 3, 1);
    printf("------------------------------------------------------------\n");

    return 0;
}
