#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"

/* Define the system of ODEs: example 2.2 (Chapter 1, pag.64) - Solution of 2-body-priblem
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.) */

typedef struct
{
    real m1;
    real m2;
    real G;

} twoBodyParams;

void twoBodyFun(const real t, const real *yy, const void *params, real *ff)
{
    // Parameters:
    assert(t < 10000000000.0); // Useless, just to compile
    twoBodyParams *Params = (twoBodyParams *)params;
    real m1 = Params->m1;
    real m2 = Params->m2;
    real G = Params->G;

    // Psotion:
    real RR1_[3] = {yy[0], yy[1], yy[2]};
    real RR2_[3] = {yy[3], yy[4], yy[5]};

    // Velocity:
    real VV1_[3] = {yy[6], yy[7], yy[8]};
    real VV2_[3] = {yy[9], yy[10], yy[11]};

    // Acceleration:
    real AA1_[3];
    real AA2_[3];

    // Compute Acceleration:
    real RR_diff[3];
    _dynorb_rvcopy(3, RR2_, RR_diff);
    _dynorb_raxpy(3, -1.0, RR1_, RR_diff);
    real r = _dynorb_rnrm2(3, RR_diff);
    _dynorb_rscal(3, (1 / (r * r * r)), RR_diff); // (R2-R1)/r^3
    _dynorb_rvcopy(3, RR_diff, AA1_);             // Compute AA1_ and AA2_
    _dynorb_rvcopy(3, RR_diff, AA2_);             // Compute AA1_ and AA2_
    _dynorb_rscal(3, (G * m2), AA1_);             // Compute AA1_ and AA2_
    _dynorb_rscal(3, (-1.0 * G * m1), AA2_);      // Compute AA1_ and AA2_

    // Update State Derivatives:
    _dynorb_rvcopy(3, VV1_, &ff[0]);
    _dynorb_rvcopy(3, VV2_, &ff[3]);
    _dynorb_rvcopy(3, AA1_, &ff[6]);
    _dynorb_rvcopy(3, AA2_, &ff[9]);
}

int main(void)
{
    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */
    // Step 1: Set up of ode-System and Solver Configuration:

    // Declare Structures:
    _dynorb_odeSys twoBodySys;
    _dynorb_solverConf solverConf;

    // Structure Fields:
    twoBodyParams p =
        {
            .m1 = 1.0e26,
            .m2 = 1.0e26,
            .G = 6.67259e-20,
        };

    real RR1_0[3] = {0, 0, 0};
    real RR2_0[3] = {3000, 0, 0};
    real VV1_0[3] = {10, 20, 30};
    real VV2_0[3] = {0, 40, 0};

    // Initial State:
    real yy0[12];
    _dynorb_rvcopy(3, RR1_0, &yy0[0]);
    _dynorb_rvcopy(3, RR2_0, &yy0[3]);
    _dynorb_rvcopy(3, VV1_0, &yy0[6]);
    _dynorb_rvcopy(3, VV2_0, &yy0[9]);

    // Solver configuration:
    int sys_size = 12; // Size of twoBodySys
    real t0 = 0.0;     // Initial Time [s]
    real t1 = 480;     // Final Time [s] // 8.0 * 60.0;
    real h = 0.1;      // Time step [s]

    // Confugure (dynamically)
    _dynorb_configure_static(&twoBodySys, &solverConf,
                             &twoBodyFun, &p, yy0,
                             sys_size, t0, t1, h);

    real tt[solverConf.n_steps];
    real YY_t[solverConf.n_steps * twoBodySys.sys_size];
    twoBodySys.tt = tt;
    twoBodySys.YY_t = YY_t;

    // Step 3: Perform Integration using custom RK4 method:
    _dynorb_rrk4(&twoBodySys, &solverConf);
    printf("\nDONE INTEGRATING\n");

    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    // Step 4: Open file and Loop over time steps save solution
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

    for (int i = 0; i < solverConf.n_steps; i++)
    {

        // [Time, RR1, RR2, VV1, VV2]
        fprintf(outfile1, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f \n",
                twoBodySys.tt[i],
                twoBodySys.YY_t[i * twoBodySys.sys_size], twoBodySys.YY_t[i * twoBodySys.sys_size + 1], twoBodySys.YY_t[i * twoBodySys.sys_size + 2],
                twoBodySys.YY_t[i * twoBodySys.sys_size + 3], twoBodySys.YY_t[i * twoBodySys.sys_size + 4], twoBodySys.YY_t[i * twoBodySys.sys_size + 5]);

        // Center of Mass:
        _dynorb_rvcopy(3, &twoBodySys.YY_t[i * twoBodySys.sys_size], RRG);          // Copy position of m1 to RRG
        _dynorb_rscal(3, p.m1, RRG);                                                // Scale by m1
        _dynorb_raxpy(3, p.m2, &twoBodySys.YY_t[i * twoBodySys.sys_size + 3], RRG); // Add m2 * RR2 to the CoM
        _dynorb_rscal(3, 1.0 / (p.m1 + p.m2), RRG);                                 // Divide by total mass (m1 + m2)

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

    // Step 5: Free Memory and Close the files:
    // _dynorb_free(&twoBodySys);
    fclose(outfile1);
    fclose(outfile2);
    fclose(outfile3);
    fclose(outfile4);
    // printf("\nDONE FREEING MEMORY AND LOGGING DATA\n");
    printf("\nDONE LOGGING DATA\n");

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command1 =
        "set terminal qt enhanced\n"
        "set title 'Example 2 Chapter 02: 3D Trajectories - Ansolute - Inertial Frame [km]'\n"
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
        "set title 'Example 2 Chapter 02: 3D Trajectories - Relative to body 1 - Inertial Frame [km]'\n"
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
        "set title 'Example 2 Chapter 02: 3D Trajectories - Relative to body 1 - Inertial Frame [km]'\n"
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

/**
 * @brief Solves the inertial two-body problem in three dimensions numerically using the RKF4(5) method.
 *
 * This function uses the RKF4(5) method to solve the gravitational two-body problem in 3D.
 * It computes the motion of two bodies under their mutual gravitational attraction.
 *
 * @param G Universal gravitational constant (km^3/kg/s^2).
 * @param m1 Mass of the first body (kg).
 * @param m2 Mass of the second body (kg).
 * @param m Total mass (kg).
 * @param t0 Initial time (s).
 * @param tf Final time (s).
 * @param R1_0 3x1 column vector containing the components of the initial position (km) of m1.
 * @param V1_0 3x1 column vector containing the components of the initial velocity (km/s) of m1.
 * @param R2_0 3x1 column vector containing the components of the initial position (km) of m2.
 * @param V2_0 3x1 column vector containing the components of the initial velocity (km/s) of m2.
 * @param y0 12x1 column vector containing the initial values of the state vectors of the two bodies: [R1_0; R2_0; V1_0; V2_0].
 * @param t Column vector of the times at which the solution is found.
 * @param X1 Column vector containing the X coordinates (km) of m1 at the times in t.
 * @param Y1 Column vector containing the Y coordinates (km) of m1 at the times in t.
 * @param Z1 Column vector containing the Z coordinates (km) of m1 at the times in t.
 * @param X2 Column vector containing the X coordinates (km) of m2 at the times in t.
 * @param Y2 Column vector containing the Y coordinates (km) of m2 at the times in t.
 * @param Z2 Column vector containing the Z coordinates (km) of m2 at the times in t.
 * @param VX1 Column vector containing the X component of the velocity (km/s) of m1 at the times in t.
 * @param VY1 Column vector containing the Y component of the velocity (km/s) of m1 at the times in t.
 * @param VZ1 Column vector containing the Z component of the velocity (km/s) of m1 at the times in t.
 * @param VX2 Column vector containing the X component of the velocity (km/s) of m2 at the times in t.
 * @param VY2 Column vector containing the Y component of the velocity (km/s) of m2 at the times in t.
 * @param VZ2 Column vector containing the Z component of the velocity (km/s) of m2 at the times in t.
 * @param y A matrix whose 12 columns are respectively [X1, Y1, Z1; X2, Y2, Z2; VX1, VY1, VZ1; VX2, VY2, VZ2].
 * @param XG Column vector containing the X coordinates (km) of the center of mass at the times in t.
 * @param YG Column vector containing the Y coordinates (km) of the center of mass at the times in t.
 * @param ZG Column vector containing the Z coordinates (km) of the center of mass at the times in t.
 *
 * @note This function requires `rk4`.
 */
