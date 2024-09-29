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
    printf("Final Position Vector (Semi-analytical) [km]:\nrr = \n");
    _dynorb_rmprint(&yy[0], 3, 1);
    printf("\n");
    printf("Final Velocity Vector (Semi-analytical) [km/s]:\nvv = \n");
    _dynorb_rmprint(&yy[3], 3, 1);
    printf("------------------------------------------------------------\n");

    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 1: System and Solver Parameters, Initial Conditions
     * ---------------------------------------------------------- */

    // Two Body System Parameters:
    _dynorb_twoBodyRelParams twoBodyParams;
    twoBodyParams.mu = mu_E; // [km^3/s^2]

    // Solver Parameters:
    int sys_size = 6;  // Size of twoBodySys
    real t0 = 0.0;     // Initial Time [s]
    real t1 = t0 + Dt; // Final Time [s] (4 hours)
    real h = 1;        // Time step [s]

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
    _dynorb_rrk4(&twoBodySys, &solverConf);

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Final Position Vector (Numerical Integration) [km]:\nrr = \n");
    _dynorb_rmprint(&twoBodySys.YY_t[(solverConf.n_steps - 1) * twoBodySys.sys_size], 3, 1);
    printf("\n");
    printf("Final Velocity Vector (Numerical Integration) [km/s]:\nvv = \n");
    _dynorb_rmprint(&twoBodySys.YY_t[(solverConf.n_steps - 1) * twoBodySys.sys_size + 3], 3, 1);
    printf("------------------------------------------------------------\n");

    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 4: Open files and log data:
     * ---------------------------------------------------------- */
    FILE *outfile1 = fopen("./data/ex_03_07_a.txt", "w");
    if (outfile1 == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    FILE *outfile2 = fopen("./data/ex_03_07_b.txt", "w");
    if (outfile2 == NULL)
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

    // Initial Time and State:
    fprintf(outfile2, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
            t0, yy0[0], yy0[1], yy0[2], yy0[3], yy0[4], yy0[5]);
    // Final Time and State:
    fprintf(outfile2, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
            t1, yy[0], yy[1], yy[2], yy[3], yy[4], yy[5]);

    /* ----------------------------------------------------------
     * STEP 5: Memory Deallocation (dynamic) and file closure:
     * ---------------------------------------------------------- */
    _dynorb_free(&twoBodySys);
    fclose(outfile1);
    fclose(outfile2);
    printf("\nDONE FREEING MEMORY AND LOGGING DATA\n");

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 6: Plot data:
     * ---------------------------------------------------------- */
    const char *plot_command1 =
        "set terminal qt enhanced\n"
        "set title 'Example 7 Chapter 03: 3D Trajectories - Relative motion - ECI Coordinates'\n"
        "set xlabel 'X [km]'\n"
        "set ylabel 'Y [km]'\n"
        "set zlabel 'Z [km]'\n"
        "set size ratio -1\n"
        "set grid\n"
        "set view 0, 0\n"
        // "set key right top\n" // Add this line to move the legend
        "set key at screen 0.999,0.999\n" // Example for finer adjustment
        "splot './data/ex_03_07_a.txt' using 2:3:4 with points pt 7 ps 0.3 lc rgb 'black' title 'orbit', \\\n"
        "      './data/ex_03_07_a.txt' using 5:6:7 with points pt 7 ps 10 lc rgb 'blue' notitle, \\\n"
        "      './data/ex_03_07_a.txt' every ::0::0 using 2:3:4 with points pt 7 ps 1 lc rgb 'green' title 'start', \\\n"
        "      './data/ex_03_07_a.txt' every ::3599::3599 using 2:3:4 with points pt 7 ps 1 lc rgb 'red' title 'end', \\\n"
        "      './data/ex_03_07_b.txt' every ::0::0 using 2:3:4 with points pt 1 ps 5 lc rgb 'green' title 'kepStart', \\\n"
        "      './data/ex_03_07_b.txt' every ::1::1 using 2:3:4 with points pt 1 ps 5 lc rgb 'red' title 'kepEnd'\n";

    FILE *gnuplot1 = popen("gnuplot -persistent", "w");
    if (gnuplot1)
    {
        fprintf(gnuplot1, "%s", plot_command1);
        pclose(gnuplot1);
    }

    return 0;
}
