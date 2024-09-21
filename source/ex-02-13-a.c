#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"
// Time integration
#include <sys/time.h>

/*  Example 2.13 (Chapter 2, pag.117): Compute the state vector [rr,vv] from the initial
    state vector [rr0,vv0] and the change in true anomaly Dth.
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

int main(void)
{
    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 1: Solution
     * ---------------------------------------------------------- */

    // Earth Gravitational Parameter:
    real mu_E = _dynorb_MU_E; // [km^3/s^2]

    // Initial Conditions:
    real rr_0[3] = {8182.4, -6865.9, 0.0}; //[km]
    real vv_0[3] = {0.47572, 8.8116, 0};   //[km/s]
    real yy0[6];                           // Initial State
    _dynorb_rvcopy(6, rr_0, &yy0[0]);
    _dynorb_rvcopy(6, vv_0, &yy0[3]);

    // Difference in True Anomaly:
    real Dth_rad = (120.0 / 180.0) * _dynorb_PI;

    // Next State:
    real yy[6];
    _dynorb_yy_From_yy0_Dth(mu_E, Dth_rad, yy0, yy);

    // Lagrange Functions:
    real fgdfdg[4];
    _dynorb_LagrangeFunctionsFrom_yy0_Dth(mu_E, Dth_rad, yy0, fgdfdg);

    // Print:
    printf("\n");
    printf("f = %f\n", fgdfdg[0]);
    printf("g = %f\n", fgdfdg[1]);
    printf("f_dot = %f\n", fgdfdg[2]);
    printf("g_dot = %f\n", fgdfdg[3]);
    printf("\n");
    printf("rr0 [km]:\n");
    _dynorb_rmprint(yy, 3, 1);
    printf("\n");
    printf("vv0 [km/s]:\n");
    _dynorb_rmprint(&yy[3], 3, 1);

    return 0;
}
