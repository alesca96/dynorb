#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

/*  Example 3.5 (Chapter 3, pag.170):This program uses Algorithm 3.2 and the data
    of Example 3.5 to solve Kepler’s equation.
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

int main(void)
{
    /* DATA: */
    real rp = (real)_dynorb_R_E + 300.0; // Perigee Radius [km]
    real vp = 15;                        // Velocity at Perigee [km/s]
    real th = _dynorb_deg2rad(100.0);    // True Anomaly [rad]
    real mu_E = (real)_dynorb_MU_E;      // Earth grav. parameter [km^3 /s^2]

    /* EXERCISE 3.5: */
    real h = rp * vp;                                                           // Specific Angular Momentum [km^2/s]
    real e = ((h * h) / (mu_E * rp)) - 1.0;                                     // Eccentricity [-]
    real th_inf = acos(-1.0 / e);                                               // True Anomaly of the asymptote [rad]
    real r_th = ((h * h) / mu_E) * (1 / (1 + e * cos(th)));                     // Radius at th = 100 deg [km]
    real F = 2.0 * atanh(sqrt((e - 1) / (e + 1)) * tan(th / 2.0));              // Hyperbolic Eccentric Anomaly [rad]
    real Mh = e * sinh(F) - F;                                                  // Mean Anomaly (hyperbola) [rad]
    real tp = (pow(h, 3.0) / (mu_E * mu_E)) * (1 / pow((e * e) - 1, 1.5)) * Mh; // Time since Perigee [s]

    // Print Info:
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.5: part (a)\n");
    printf("-----------------------\n\n");
    printf("True Anomaly th = %g [deg]\n", _dynorb_rad2deg(th));
    printf("Angular Momentum h = %gm [km^2/s]\n", h);
    printf("Eccentricity e = %g [-]\n", e);
    printf("Radius at th = %g [deg] r(th) = %g [km]\n", _dynorb_rad2deg(th), r_th);
    printf("True Anomaly of Asymptote th_inf = %g [deg]\n", _dynorb_rad2deg(th_inf));
    printf("Hyperbolic Eccentric Anomaly (from Kepler) F = %g [rad]\n", F);
    printf("Mean Anomaly Mh= %g [rad]\n", (Mh));
    printf("Time since Perigee tp = %g [s]\n", tp);
    printf("------------------------------------------------------------\n");

    // Part b:
    real tp_ = tp + (3.0 * 60.0 * 60.0);                              // 3 hours after previous true anomaly (th = 100 deg) [s]
    Mh = ((mu_E * mu_E) / pow(h, 3.0)) * pow((e * e) - 1, 1.5) * tp_; // [rad]
    F = _dynorb_rkeplerH(e, Mh, 1.e-6);                               //  Hyperbolic Eccentric Anomaly [rad]
    real th_ = 2.0 * atan(sqrt((e + 1) / (e - 1)) * tanh(F / 2.0));   //  New True Anomaly [rad]
    real r_th_ = ((h * h) / mu_E) * (1 / (1 + e * cos(th_)));         // [km]
    real v_tan = (h / r_th_);                                         // Tangential Velocity [km/s]
    real v_rad = (mu_E / h) * e * sin(th_);                           // Radial Velocity [km/s]
    real v_th_ = sqrt((v_tan * v_tan) + (v_rad * v_rad));             // Orbital Speed [km/s]
    real v_inf = (mu_E / h) * e * sin(th_inf);                        // Hyperbolic Eccess Speed [km/s]

    // Print Info
    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Example 3.5: part (b)\n");
    printf("-----------------------\n\n");
    printf("New Time since Perigee tp = %g [s]\n", tp_);
    printf("Mean Anomaly Mh= %g [rad]\n", (Mh));
    printf("Hyperbolic Eccentric Anomaly (from Kepler) F = %g [rad]\n", F);
    printf("New True Anomaly th = %g [deg]\n", _dynorb_rad2deg(th_));
    printf("Radius at th = %g [deg] r(th) = %g [km]\n", _dynorb_rad2deg(th_), r_th_);
    printf("Speed at th = %g [deg] v(th) = %g [km/s]\n", _dynorb_rad2deg(th_), v_th_);
    printf("Hyperbolic Ecces Speed at th_inf = %g [deg] v(th) = %g [km/s]\n", _dynorb_rad2deg(th_inf), v_inf);
    printf("------------------------------------------------------------\n");

    return 0;
}
