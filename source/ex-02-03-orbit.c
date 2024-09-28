#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"
// Time integration
#include <sys/time.h>

/* Define the system of ODEs: example 2.2 (Chapter 1, pag.64) - Solution of 2-body-priblem
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.) */

int main(void)
{
    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 1: System and Solver Parameters, Initial Conditions
     * ---------------------------------------------------------- */

    // Two Body System Parameters:
    real m_E = 5.974e24; // [kg]
    real G = 6.6742e-20; // [km^3/kg/s^2]
    _dynorb_twoBodyRelParams twoBodyParams;
    twoBodyParams.m = 1000.0;                       // [kg]
    twoBodyParams.mu = G * (twoBodyParams.m + m_E); // [km^3/s^2]

    // Solver Parameters:
    int sys_size = 6;      // Size of twoBodySys
    real t0 = 0.0;         // Initial Time [s]
    real t1 = 4 * 60 * 60; // Final Time [s] (4 hours)
    real h = 1;            // Time step [s]

    // Initial Conditions:
    real rr_0[3] = {8000, 0, 6000}; //[km]
    real vv_0[3] = {0, 7, 0};       //[km/s]
    real yy0[6];                    // Initial State
    _dynorb_rvcopy(6, rr_0, &yy0[0]);
    _dynorb_rvcopy(6, vv_0, &yy0[3]);

    /* ----------------------------------------------------------
     * STEP 2: Confugure ode-system and solver:
     * ---------------------------------------------------------- */

    // Structures declaration:
    _dynorb_odeSys twoBodySys;
    _dynorb_solverConf solverConf;
    // Configuration
    _dynorb_configure_dynamic(&twoBodySys, &solverConf,
                              &_dynorb_twoBodyRelFun, &twoBodyParams, yy0,
                              sys_size, t0, t1, h);

    /* ----------------------------------------------------------
     * STEP 3: Numerical Integration:
     * ---------------------------------------------------------- */
    struct timeval start, end;
    long seconds, useconds;
    double time_taken;
    gettimeofday(&start, NULL);
    _dynorb_rrk4(&twoBodySys, &solverConf);
    gettimeofday(&end, NULL);
    // Calculate the time difference in seconds and microseconds
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    time_taken = (seconds * 1.e6) + useconds;
    printf("\nDONE INTEGRATING\n");
    printf("Elapsed time: %f micro-seconds\n", time_taken);

    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 4: Open files and log data:
     * ---------------------------------------------------------- */
    FILE *outfile1 = fopen("./data/ex_02_03_a.txt", "w");
    if (outfile1 == NULL)
    {
        perror("Error opening file");
        return 1;
    }
    // Loop over data to log it:
    for (int i = 0; i < solverConf.n_steps; i++)
    {
        // [Time, rr1, rr2, vv, vv]
        fprintf(outfile1, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
                twoBodySys.tt[i],
                twoBodySys.YY_t[i * twoBodySys.sys_size], twoBodySys.YY_t[i * twoBodySys.sys_size + 1], twoBodySys.YY_t[i * twoBodySys.sys_size + 2],
                twoBodySys.YY_t[i * twoBodySys.sys_size + 3], twoBodySys.YY_t[i * twoBodySys.sys_size + 4], twoBodySys.YY_t[i * twoBodySys.sys_size + 5]);
    }

    /* ----------------------------------------------------------
     * STEP 5: Memory Deallocation (dynamic) and file closure:
     * ---------------------------------------------------------- */
    _dynorb_free(&twoBodySys);
    fclose(outfile1);

    printf("\nDONE FREEING MEMORY AND LOGGING DATA\n");
    // printf("\nDONE LOGGING DATA\n");

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 6: Plot data:
     * ---------------------------------------------------------- */
    const char *plot_command1 =
        "set terminal qt enhanced\n"
        "set title 'Example 2 Chapter 02: 3D Trajectories - Relative motion - ECI Coordinates'\n"
        "set xlabel 'X [km]'\n"
        "set ylabel 'Y [km]'\n"
        "set zlabel 'Z [km]'\n"
        "set size ratio -1\n"
        "set grid\n"
        "set view 60, 330\n"
        "splot './data/ex_02_03_a.txt' using 2:3:4 with points pt 7 ps 0.3 lc rgb 'black' title 'orbit', \\\n"
        "      './data/ex_02_03_a.txt' using 5:6:7 with points pt 7 ps 10 lc rgb 'blue' notitle, \\\n"
        "      './data/ex_02_03_a.txt' every ::0::0 using 2:3:4 with points pt 7 ps 1 lc rgb 'green' title 'start', \\\n"
        "      './data/ex_02_03_a.txt' every ::14390::14390 using 2:3:4 with points pt 7 ps 1 lc rgb 'red' title 'end'\n";

    FILE *gnuplot1 = popen("gnuplot -persistent", "w");
    if (gnuplot1)
    {
        fprintf(gnuplot1, "%s", plot_command1);
        pclose(gnuplot1);
    }

    return 0;
}
