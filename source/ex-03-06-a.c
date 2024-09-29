#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

/*  Example 3.5 (Chapter 3, pag.170):This program uses Algorithm 3.2 and the data
    of Example 3.5 to solve Kepler’s equation.
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

int main(void)
{
    // DATA:
    real mu_E = (real)_dynorb_MU_E;   // Earth grav. parameter [km^3 /s^2]
    real th0 = _dynorb_deg2rad(30.0); // True Anomaly [rad]
    real r0 = 10000;                  // Initial Radial Distance [km]
    real vr0 = 3.0752;                // Initial Radial Velocity [km/s]
    real e = 1.4682;                  // Eccentricity [-]
    real Dt = 1.0 * 60 * 60;          // Time Passed from th0 [s]

    // EXERCISE 3.1:
    real h = sqrt(mu_E * r0 * (1 + e * cos(th0)));                             // Specific Angular Momentum [km^2/s]
    real th_inf = acos(-1.0 / e);                                              // True Anomaly of the asymptote [rad]
    real r_th0 = ((h * h) / mu_E) * (1 / (1 + e * cos(th0)));                  // Radius at th0 = 100 deg [km]
    real F0 = 2.0 * atanh(sqrt((e - 1) / (e + 1)) * tan(th0 / 2.0));           // Hyperbolic Eccentric Anomaly [rad]
    real a = ((h * h) / mu_E) * (1 / (1 - (e * e)));                           // Hyperbola Semi-major Axis [km]
    real alpha = 1 / a;                                                        // Alpha
    real chi = _dynorb_rkeplerUniversal(mu_E, Dt, r0, vr0, alpha, 1.e-6, 100); // Universal Variable [km^0.5]
    real F = F0 + chi / sqrt(-a);                                              // New  Hyperbolic Eccentric Anomaly [rad]
    real th = 2.0 * atan(sqrt((e + 1) / (e - 1)) * tanh(F / 2.0));             //  New True Anomaly [rad]

    // Print Info:
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.6: \n");
    printf("-----------------------\n\n");
    printf("True Anomaly th0 = %g [deg]\n", _dynorb_rad2deg(th0));
    printf("Angular Momentum h = %gm [km^2/s]\n", h);
    printf("Eccentricity e = %g [-]\n", e);
    printf("Radius at th0 = %g [deg]: r(th0) = %g [km]\n", _dynorb_rad2deg(th0), r_th0);
    printf("Radial Velocity at th0 = %g [deg]: vr(th0) = %g [km/s]\n", _dynorb_rad2deg(th0), vr0);
    printf("True Anomaly of Asymptote th_inf = %g [deg]\n", _dynorb_rad2deg(th_inf));
    printf("Semi-major Axis a = %g [km]\n", a);
    printf("Hyperbolic Eccentric Anomaly (from Kepler) F0 = %g [rad]\n", F0);
    printf("Universal Variable Chi = %g [km^0.5]\n", chi);
    printf(" New Hyperbolic Eccentric Anomaly F = %g [rad]\n", F);
    printf("New True Anomaly th = %g [deg]\n", _dynorb_rad2deg(th));
    printf("------------------------------------------------------------\n");

    return 0;
}

/*
int main(void)
{
    // DATA:
    real mu_E = (real)_dynorb_MU_E; // Earth grav. parameter [km^3 /s^2]
    real r0 = 10000;                // Initial Radial Distance [km]
    real vr0 = 3.0752;              // Initial Radial Velocity [km/s]
    real a = -19655;                // Hyperbola Semi-major Axis [km]
    real Dt = 1.0 * 60 * 60;        // Time Passed from th0 [s]

    // EXERCISE 3.1:
    real alpha = 1 / a;                                                        // Alpha
    real chi = _dynorb_rkeplerUniversal(mu_E, Dt, r0, vr0, alpha, 1.e-6, 100); // Universal Variable [km^0.5]

    // Print Info:
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.6: \n");
    printf("-----------------------\n\n");
    printf("Initial Radial Coordinate r0 = %g [km]\n", r0);
    printf("Initial Radial Velocity r0 = %g [km/s]\n", vr0);
    printf("Elapsed Time Dt = %g [s]\n", Dt);
    printf("Semi-major Axis a = %g [km]\n", a);
    printf("Universal Variable Chi = %g [km^0.5]\n", chi);
    printf("------------------------------------------------------------\n");

    return 0;
}
*/
