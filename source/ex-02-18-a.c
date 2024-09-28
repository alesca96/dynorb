#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

/*  Example 2.18 (Chapter 2, pag.137):This program uses the RK4 method to solve the
    earth-moon restricted three-body problem (Equations 2.192a and 2.192b) for the
    trajectory of a spacecraft having the initial specified conditions.
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

int main(void)
{
    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 1: System and Solver Parameters, Initial Conditions
     * ---------------------------------------------------------- */

    /* CONSTANTS and PARAMETERS: */
    real sec_in_day = 24 * 60 * 60; // Seconds in a day [s]
    real G = 6.6742e-20;            // Unversal grav. const. [km^3/(kg s^2)]
    // real rmoon = 1737;                     // Moon radius [km]
    real rearth = 6378;                    // Earth radis [km]
    real r12 = 384400;                     // Earth-Moon distance [km]
    real m1 = 5974e21;                     // Earth mass [kg]
    real m2 = 7348e19;                     // Moon mass [kg]
    real M = m1 + m2;                      // Total mass [kg]
    real pi_1 = m1 / M;                    // Earth mass ratio [-]
    real pi_2 = m2 / M;                    // Moon mass ratio [-]
    real mu1 = G * m1;                     // Earth grav. parameter [km^3 /s^2]
    real mu2 = G * m2;                     // Moon grav. parameter [km^3 /s^2]
    real mu = mu1 + mu2;                   // Earth-Moon system grav. parameter [km^3 /s^2]
    real W = sqrt(mu / (r12 * r12 * r12)); // Angular velocity Moon around Earth [rad/s]
    real x1 = -1.0 * pi_2 * r12;           // x-coordinates of the Earth relative to the Earth-Moon barycenter [km]
    real x2 = 1.0 * pi_1 * r12;            // x-coordinates of the Earth relative to the Earth-Moon barycenter [km]

    /* INITIAL CONDITIONS and INPUTS: */
    real d0 = 200.0;                                                 // Initial altitude s/c [km]
    real phi = _dynorb_deg2rad(-90);                                 // Polar Azimuth coordinate of s/c (+ counterclockwise from Earth-Monn line) [rad]
    real v0 = 10.9148;                                               // initial speed of spacecraft relative to rotating Earth-Moon system [km/s]
    real gamma = _dynorb_deg2rad(20);                                // Initial flight path angle [rad]
    real t0 = 0;                                                     // Initial time [s]
    real tf = 3.16689 * sec_in_day;                                  // Final time [s]
    real r0 = rearth + d0;                                           // Initial radial dist. of s/c from Earth [km]
    real x0 = r0 * cos(phi) + x1;                                    // x coordinate of s/c in rotating Earth-Moon sys [km]
    real y0 = r0 * sin(phi);                                         // y coordinate of s/c in rotating Earth-Moon sys [km]
    real vx0 = v0 * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi)); // x coordinate of s/c velocity in rotating Earth-Moon sys [km/s]
    real vy0 = v0 * (sin(gamma) * sin(phi) - cos(gamma) * cos(phi)); // y coordinate of s/c velocity in rotating Earth-Moon sys [km/s]

    // Parameter Structure:
    _dynorb_threeBodyRestrictParams p;
    p.mu1 = mu1;
    p.mu2 = mu2;
    p.pi_1 = pi_1;
    p.pi_2 = pi_2;
    p.W = W;
    p.x1 = x1;
    p.x2 = x2;

    // Solver:
    int sys_size = 4; // Size of threeBodyRestrictSys
    real h = 60;      // Time step [s]

    // Initial State:
    real yy0[4] = {x0, y0, vx0, vy0};

    /* ----------------------------------------------------------
     * STEP 2: Confugure ode-system and solver:
     * ---------------------------------------------------------- */
    _dynorb_odeSys threeBodyRestrictSys;
    _dynorb_solverConf solverConf;
    _dynorb_configure_dynamic(&threeBodyRestrictSys, &solverConf,
                              &_dynorb_threeBodyRestrictFun, &p, yy0,
                              sys_size, t0, tf, h);

    /* ----------------------------------------------------------
     * STEP 3: Numerical Integration:
     * ---------------------------------------------------------- */
    _dynorb_rrk4(&threeBodyRestrictSys, &solverConf);
    printf("\nDONE INTEGRATING\n");

    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 4: Open files and log data:
     * ---------------------------------------------------------- */
    FILE *outfile = fopen("./data/ex_02_18_a.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }
    for (int i = 0; i < solverConf.n_steps; i++)
    {
        // Here Time is also recaled to minutes
        fprintf(outfile, "%.20f %.20f %.20f\n", threeBodyRestrictSys.tt[i],
                threeBodyRestrictSys.YY_t[i * threeBodyRestrictSys.sys_size],
                threeBodyRestrictSys.YY_t[i * threeBodyRestrictSys.sys_size + 1]);
    }

    /* ----------------------------------------------------------
     * STEP 5: Memory Deallocation (dynamic) and file closure:
     * ---------------------------------------------------------- */
    _dynorb_free(&threeBodyRestrictSys);
    fclose(outfile);
    printf("\nDONE FREEING MEMORY AND LOGGING DATA\n");

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 6: Plot data:
     * ---------------------------------------------------------- */
    const char *plot_command1 =
        "set terminal qt enhanced\n"
        "set title 'Example 18 Chapter 02: Mass under GravityRestricted 3-body problem - Position [km]'\n"
        "set xlabel 'x(t) [km]'\n"
        "set ylabel 'y(t) [km]'\n"
        "set grid\n"
        "plot './data/ex_02_18_a.txt' using 2:3 with lines lc rgb 'red' title 'r(t)'\n";

    FILE *gnuplot1 = popen("gnuplot -persistent", "w");
    if (gnuplot1)
    {
        fprintf(gnuplot1, "%s", plot_command1);
        pclose(gnuplot1);
    }

    return 0;
}

/*
tt: time ticks vector
ff: matrix with
 column 1: solution for x at the times in t
 column 2: solution for y at the times in t
 column 3: solution for vx at the times in t
 column 4: solution for vy at the times in t
 xf,yf- x and y coordinates of spacecraft in rotating earth-moon
 system at tf
 vxf, vyf- x and y components of spacecraft velocity relative to
 rotating earth-moon system at tf
 df- distance from surface of the moon at tf
 vf- relative speed at tf
 */