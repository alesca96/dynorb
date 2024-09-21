#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

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
    _dynorb_twoBodyAbsParams twoBodyParams;
    twoBodyParams.m1 = 1.0e26;
    twoBodyParams.m2 = 1.0e26;
    twoBodyParams.G = 6.67259e-20;

    // Solver Parameters:
    int sys_size = 12; // Size of twoBodySys
    real t0 = 0.0;     // Initial Time [s]
    real t1 = 480;     // Final Time [s] // 8.0 * 60.0;
    real h = 0.1;      // Time step [s]

    // Initial Conditions:
    real RR1_0[3] = {0, 0, 0};
    real RR2_0[3] = {3000, 0, 0};
    real VV1_0[3] = {10, 20, 30};
    real VV2_0[3] = {0, 40, 0};
    real yy0[12]; // Initial State:
    _dynorb_rvcopy(3, RR1_0, &yy0[0]);
    _dynorb_rvcopy(3, RR2_0, &yy0[3]);
    _dynorb_rvcopy(3, VV1_0, &yy0[6]);
    _dynorb_rvcopy(3, VV2_0, &yy0[9]);

    /* ----------------------------------------------------------
     * STEP 2: Confugure ode-system and solver:
     * ---------------------------------------------------------- */

    // Structures declaration:
    _dynorb_odeSys twoBodySys;
    _dynorb_solverConf solverConf;
    // Configuration
    _dynorb_configure_dynamic(&twoBodySys, &solverConf,
                              &_dynorb_twoBodyAbsFun, &twoBodyParams, yy0,
                              sys_size, t0, t1, h);
    // // Static memory allocation:
    // real tt[solverConf.n_steps];
    // real YY_t[solverConf.n_steps * twoBodySys.sys_size];
    // twoBodySys.tt = tt;
    // twoBodySys.YY_t = YY_t;

    /* ----------------------------------------------------------
     * STEP 3: Numerical Integration:
     * ---------------------------------------------------------- */
    _dynorb_rrk4(&twoBodySys, &solverConf);
    printf("\nDONE INTEGRATING\n");

    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    /* ----------------------------------------------------------
     * STEP 4: Open files and log data:
     * ---------------------------------------------------------- */
    FILE *outfile1 = fopen("./data/ex_02_02_a.txt", "w");
    FILE *outfile2 = fopen("./data/ex_02_02_a1_CoM.txt", "w");
    FILE *outfile3 = fopen("./data/ex_02_02_a2_RelToBody1.txt", "w");
    FILE *outfile4 = fopen("./data/ex_02_02_a3_RelToCoM.txt", "w");
    if (outfile1 == NULL || outfile2 == NULL || outfile3 == NULL || outfile4 == NULL)
    {
        perror("Error opening file");
        return 1;
    }
    // Center of Mass Position:
    real RRG[3];
    // Relative Positions:
    real RR1_rel[3];
    real RR2_rel[3];
    real RRG_rel[3];
    // Loop over data to log it:
    for (int i = 0; i < solverConf.n_steps; i++)
    {

        // [Time, RR1, RR2, VV1, VV2]
        fprintf(outfile1, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
                twoBodySys.tt[i],
                twoBodySys.YY_t[i * twoBodySys.sys_size], twoBodySys.YY_t[i * twoBodySys.sys_size + 1], twoBodySys.YY_t[i * twoBodySys.sys_size + 2],
                twoBodySys.YY_t[i * twoBodySys.sys_size + 3], twoBodySys.YY_t[i * twoBodySys.sys_size + 4], twoBodySys.YY_t[i * twoBodySys.sys_size + 5]);

        // Center of Mass:
        _dynorb_rvcopy(3, &twoBodySys.YY_t[i * twoBodySys.sys_size], RRG);                      // Copy position of m1 to RRG
        _dynorb_rscal(3, twoBodyParams.m1, RRG);                                                // Scale by m1
        _dynorb_raxpy(3, twoBodyParams.m2, &twoBodySys.YY_t[i * twoBodySys.sys_size + 3], RRG); // Add m2 * RR2 to the CoM
        _dynorb_rscal(3, 1.0 / (twoBodyParams.m1 + twoBodyParams.m2), RRG);                     // Divide by total mass (m1 + m2)

        // [Time, RRG]
        fprintf(outfile2, "%.20f %.20f %.20f %.20f \n",
                twoBodySys.tt[i],
                RRG[0], RRG[1], RRG[2]);

        // Postions relative to Mass 1:
        _dynorb_rvcopy(3, &twoBodySys.YY_t[i * twoBodySys.sys_size], RR1_rel);
        _dynorb_rvcopy(3, &twoBodySys.YY_t[i * twoBodySys.sys_size + 3], RR2_rel);
        _dynorb_rvcopy(3, RRG, RRG_rel);
        _dynorb_raxpy(3, -1.0, &twoBodySys.YY_t[i * twoBodySys.sys_size], RR1_rel);
        _dynorb_raxpy(3, -1.0, &twoBodySys.YY_t[i * twoBodySys.sys_size], RR2_rel);
        _dynorb_raxpy(3, -1.0, &twoBodySys.YY_t[i * twoBodySys.sys_size], RRG_rel);

        // [Time, RR_rel]
        fprintf(outfile3, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
                twoBodySys.tt[i],
                RR1_rel[0], RR1_rel[1], RR1_rel[2],
                RR2_rel[0], RR2_rel[1], RR2_rel[2],
                RRG_rel[0], RRG_rel[1], RRG_rel[2]);

        // Postions relative to Mass 1:
        _dynorb_rvcopy(3, &twoBodySys.YY_t[i * twoBodySys.sys_size], RR1_rel);
        _dynorb_rvcopy(3, &twoBodySys.YY_t[i * twoBodySys.sys_size + 3], RR2_rel);
        _dynorb_rvcopy(3, RRG, RRG_rel);
        _dynorb_raxpy(3, -1.0, RRG, RR1_rel);
        _dynorb_raxpy(3, -1.0, RRG, RR2_rel);
        _dynorb_raxpy(3, -1.0, RRG, RRG_rel);

        // [Time, RR_rel]
        fprintf(outfile4, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
                twoBodySys.tt[i],
                RR1_rel[0], RR1_rel[1], RR1_rel[2],
                RR2_rel[0], RR2_rel[1], RR2_rel[2],
                RRG_rel[0], RRG_rel[1], RRG_rel[2]);
    }

    /* ----------------------------------------------------------
     * STEP 5: Memory Deallocation (dynamic) and file closure:
     * ---------------------------------------------------------- */
    _dynorb_free(&twoBodySys);
    fclose(outfile1);
    fclose(outfile2);
    fclose(outfile3);
    fclose(outfile4);
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
        "set title 'Example 2 Chapter 02: 3D Trajectories - Absolute motion - Inertial Coordinates'\n"
        "set xlabel 'X [km]'\n"
        "set ylabel 'Y [km]'\n"
        "set zlabel 'Z [km]'\n"
        "set grid\n"
        "set view 60, 330\n"
        "splot './data/ex_02_02_a.txt' using 2:3:4 with points pt 7 ps 0.3 lc rgb 'red' title 'm1 trajectory', \\\n"
        "      './data/ex_02_02_a.txt' using 5:6:7 with points pt 7 ps 0.3 lc rgb 'green' title 'm2 trajectory', \\\n"
        "      './data/ex_02_02_a1_CoM.txt' using 2:3:4 with points pt 7 ps 0.3 lc rgb 'blue' title 'CoM trajectory'\n";

    const char *plot_command2 =
        "set terminal qt enhanced\n"
        "set title 'Example 2 Chapter 02: 3D Trajectories - Relative to body 1 - Inertial Coordinates'\n"
        "set xlabel 'X [km]'\n"
        "set ylabel 'Y [km]'\n"
        "set zlabel 'Z [km]'\n"
        "set grid\n"
        "set view 50, 330\n"
        "splot './data/ex_02_02_a2_RelToBody1.txt' using 2:3:4 with points pt 7 ps 0.3 lc rgb 'red' title 'm1 trajectory', \\\n"
        "      './data/ex_02_02_a2_RelToBody1.txt' using 5:6:7 with points pt 7 ps 0.3 lc rgb 'green' title 'm2 trajectory', \\\n"
        "      './data/ex_02_02_a2_RelToBody1.txt' using 8:9:10 with points pt 7 ps 0.3 lc rgb 'blue' title 'CoM trajectory'\n";

    const char *plot_command3 =
        "set terminal qt enhanced\n"
        "set title 'Example 2 Chapter 02: 3D Trajectories - Relative to CoM - Inertial Coordinates'\n"
        "set xlabel 'X [km]'\n"
        "set ylabel 'Y [km]'\n"
        "set zlabel 'Z [km]'\n"
        "set grid\n"
        "set view 50, 330\n"
        "splot './data/ex_02_02_a3_RelToCoM.txt' using 2:3:4 with points pt 7 ps 0.3 lc rgb 'red' title 'm1 trajectory', \\\n"
        "      './data/ex_02_02_a3_RelToCoM.txt' using 5:6:7 with points pt 7 ps 0.3 lc rgb 'green' title 'm2 trajectory', \\\n"
        "      './data/ex_02_02_a3_RelToCoM.txt' using 8:9:10 with points pt 7 ps 0.3 lc rgb 'blue' title 'CoM trajectory'\n";

    FILE *gnuplot1 = popen("gnuplot -persistent", "w");
    if (gnuplot1)
    {
        fprintf(gnuplot1, "%s", plot_command1);
        pclose(gnuplot1);
    }

    FILE *gnuplot2 = popen("gnuplot -persistent", "w");
    if (gnuplot2)
    {
        fprintf(gnuplot2, "%s", plot_command2);
        pclose(gnuplot2);
    }

    FILE *gnuplot3 = popen("gnuplot -persistent", "w");
    if (gnuplot3)
    {
        fprintf(gnuplot3, "%s", plot_command3);
        pclose(gnuplot3);
    }

    return 0;
}
