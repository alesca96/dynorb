#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE // USE_FLOAT // USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h" // Include your custom RK header file

/* Define the system of ODEs: example 1.20 (Chapter 1, pag.51)
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.) */

typedef struct
{
    real m;
    real mu_E;

} odeFunParams;

void odeFun(const real t, const real *yy, const void *params, real *ff)
{
    // Parameters:
    assert(t < 10000.0);
    odeFunParams *Params = (odeFunParams *)params;
    real m = Params->m;
    real mu_E = Params->mu_E;

    // Simple Harmonic Oscillator: ff = dyy/dt
    ff[0] = yy[1];
    ff[1] = -1.0 * (m * mu_E) / (yy[0] * yy[0]);
}

int main(void)
{
    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */
    // Step 1: Set up of ode-System and Solver Configuration:

    // Declare Structures:
    _dynorb_odeSys odeSys;
    _dynorb_solverConf solverConf;

    // Structure Fields:
    odeFunParams p =
        {
            .m = 1.0,
            .mu_E = 398600.0, // [Km ^3 / s^2]
        };
    real x0 = 6500.0;           // Initial Position [km]
    real x_dot0 = 7.8;          // Initial Velocity [km/s]
    real yy0[2] = {x0, x_dot0}; // Initial State
    int sys_size = 2;           // Size of odeSys
    real t0 = 0.0;              // Initial Time [s]
    real t1 = 70.0 * 60.0;      // Final Time [s]
    real h = (1 / 6.0) * 60;    // Time step [s] // 1.5 * 60; --> program crashes

    // Confugure (dynamically)
    _dynorb_configure_dynamic(&odeSys, &solverConf,
                              &odeFun, &p, yy0,
                              sys_size, t0, t1, h);

    // Step 3: Perform Integration using custom RK method:
    const real tol = 1.e-3;
    _dynorb_rrkf45(&odeSys, &solverConf, tol);

    printf("\nDONE INTEGRATING\n");
    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    // Step 4: Open file and Loop over time steps save solution
    FILE *outfile = fopen("./data/ex_01_20_a.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }
    for (int i = 0; i < solverConf.n_steps; i++)
    {
        // Here Time is also recaled to minutes
        fprintf(outfile, "%.20f %.20f %.20f %.20f\n", odeSys.tt[i], odeSys.tt[i] / 60.0, odeSys.YY_t[i * odeSys.sys_size], odeSys.YY_t[i * odeSys.sys_size + 1]);
    }
    // Step 5: Free Memory and Closing the file:
    _dynorb_free(&odeSys);
    fclose(outfile);
    printf("\nDONE FREEING MEMORY AND LOGGING DATA\n");

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command1 =
        "set terminal qt enhanced\n"
        "set title 'Example 20 Chapter 01: Mass under Gravity - Position [km]'\n"
        "set xlabel 'Time t [min]'\n"
        "set ylabel 'x(t) [km]'\n"
        "set grid\n"
        "plot './data/ex_01_20_a.txt' using 2:3 with lines lc rgb 'red' title 'x(t)'\n";

    const char *plot_command2 =
        "set terminal qt enhanced\n"
        "set title 'Example 20 Chapter 01: Mass under Gravity - Velocity [km/s]'\n"
        "set xlabel 'Time t [min]'\n"
        "set ylabel 'dx/dt (t) [km/s]'\n"
        "set grid\n"
        "plot './data/ex_01_20_a.txt' using 2:4 with lines lc rgb 'red' title 'dx/dt (t)'\n";

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

    return 0;
}
