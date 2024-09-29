#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

/*  Example 3.1 and 3.2 (Chapter 3, pag.154):This program uses Algorithm 3.1 and the data
    of Example 3.1 and 3.2 to solve Kepler’s equation.
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

int main(void)
{
    /* DATA: */
    real rp = 9600.0;                 // Perigee Radius [km]
    real ra = 21000.0;                // Apogee Radius [km]
    real th = _dynorb_deg2rad(120.0); // True Anomaly [rad]
    real mu_E = (real)_dynorb_MU_E;   // Earth grav. parameter [km^3 /s^2]

    /* EXERCISE 3.1: */
    real e = (ra - rp) / (ra + rp);                                                // Eccentricity [-]
    real h = sqrt(rp * mu_E * (1 + e));                                            // Specific Angular Momentum [km^2/s]
    real T = (2.0 * _dynorb_PI / (mu_E * mu_E)) * pow(h / sqrt(1 - (e * e)), 3.0); // Orbital Period [s]
    real E1 = 2.0 * atan(sqrt((1 - e) / (1 + e)) * tan(th / 2.0));                 // Eccentric Anomaly [rad]
    real M = E1 - e * sin(E1);                                                     // Mean anomaly (elliptic orbit) [rad]
    real tp = (M / (2.0 * _dynorb_PI)) * T;                                        // Time since Perigee [s]

    /* EXERCISE 3.2: */
    tp = 10800.0;                                            // Time since Perigee [s]
    M = 2.0 * _dynorb_PI * (tp / T);                         // Mean anomaly (elliptic orbit) [rad]
    real E = _dynorb_rkeplerE(e, M, 1.e-6);                  // Eccentric Anomaly (Kepler's eq.) [rad]
    th = 2.0 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2.0)); // Eccentric Anomaly [rad]

    // Print Info:
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.1: solving Kepler's eq. for M\n");
    printf("Perigee Radius rp = %g [km]\n", rp);
    printf("Apogee Radius ra = %g [km]\n", ra);
    printf("True Anomaly th = %g [km]\n", ra);
    printf("Eccentricity e = %g [-]\n", e);
    printf("Angular Momentum h = %g [km^2/s]\n", h);
    printf("Period T = %g [s]\n", T);
    printf("Eccentric Anomaly (from th) E1 = %g [rad]\n", (E1));
    printf("Mean Anomaly Me = %g [rad]\n", (M));
    printf("Mean Anomaly Me = %g [deg]\n", _dynorb_rad2deg(M));
    printf("Time from Perigee tp = %g [s]\n", tp);
    printf("------------------------------------------------------------\n");

    // Print Info:
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.2: solving Kepler's eq. for E\n");
    printf("Eccentricity e = %g [-]\n", e);
    printf("Mean Anomaly Me = %g [rad]\n", (M));
    printf("Mean Anomaly Me = %g [deg]\n", _dynorb_rad2deg(M));
    printf("Eccentric Anomaly (from Kepler) E = %g [rad]\n", (E));
    printf("Eccentric Anomaly (from Kepler) E = %g [deg]\n", _dynorb_rad2deg(E));
    printf("True Anomaly th = %g [rad]\n", (th));
    printf("True Anomaly th = %g [deg]\n", _dynorb_rad2deg(th));
    printf("------------------------------------------------------------\n");

    return 0;
}
