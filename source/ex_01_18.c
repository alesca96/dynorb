#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <math.h>

/* Define the system of ODEs: example 1.18 (Chapter 1, pag.45)
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.) */

typedef struct
{
    double F0;
    double m;
    double om_n;
    double zeta;
    double om;

} SimpHarmOscParams;

int SimpHarmOsc(double t, const double yy[], double ff[], void *params)
{
    // Parameters;
    SimpHarmOscParams *p = (SimpHarmOscParams *)params;
    double F0 = p->F0;
    double m = p->m;
    double om_n = p->om_n;
    double zeta = p->zeta;
    double om = p->om;

    // Simple Harmonic Oscillator: ff = df/dt
    ff[0] = yy[1];
    ff[1] = (F0 / m) * sin(om * t) - 2 * zeta * om_n * yy[1] - (om_n * om_n) * yy[0];

    return GSL_SUCCESS;
}

int main(void)
{

    /* ==========================================================
     * GSL: Use GSL to Integrate the ODE system
     * ========================================================== */

    // Step 0: Define Parameters of ODE system:
    SimpHarmOscParams p = {0};
    p.F0 = 1.0;
    p.m = 1.0;
    p.om_n = 1.0;
    p.om = 0.4 * p.om_n;
    p.zeta = 0.03;

    // Step 1: Define the ODE sys:
    gsl_odeiv2_system sys = {SimpHarmOsc, NULL, 2, &p};

    // Step 2: Initial conditions and Integration Range:
    double t = 0.0;              // Initial t
    double t1 = 110.0;           // End t
    double x0 = 0.0;             // Initial Position
    double x_dot0 = 0.0;         // Initial Velocity
    double yy[2] = {x0, x_dot0}; // Initial State
    const double h = 0.5;        // Step size
    const double abstol = 1e-6;  // Absolute Tolerance
    const double reltol = 1e-12; // Relative Tolerance

    // Step 3: Choose a Solver (RK4) and setup stepping function:
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk2, h, abstol, reltol);

    // Step 4: Open file to store results
    FILE *outfile = fopen("./data/ex_01_18.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Step 5: Perform Integration:
    while (t < t1)
    {
        double next_t = t + h; // Calculate the next time step
        if (next_t > t1)
        {
            next_t = t1; // Ensure we do not exceed the final time
        }

        // printf("Time istant before integration: t = %f\n", t);
        int status = gsl_odeiv2_driver_apply(driver, &t, next_t, yy);
        // printf("Time istant after integration: t = %f\n", t);

        if (status != GSL_SUCCESS)
        {
            printf("Error: %s\n", gsl_strerror(status));
            break;
        }
        // printf("t = %.5f, y1 = %.5f, y2 = %.5f\n", t, yy[0], yy[1]); // Print results
        fprintf(outfile, "%f %f %f\n", t, yy[0], yy[1]);
    }

    // Step 6: Free Memory and Close Data File:
    gsl_odeiv2_driver_free(driver);
    fclose(outfile);

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command =
        "set terminal qt\n"
        "set title 'Example 18 Chapter 01: Simple Harmonic Oscillator using GSL'\n"
        "set xlabel 'Time t [s]'\n"
        "set ylabel 'x(t) [m], v(t) [m/s]'\n"
        "plot './data/ex_01_18.txt' using 1:2 with lines title 'x(t)', "
        "'./data/ex_01_18.txt' using 1:3 with lines title 'v(t)'\n";

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    return 0;
}
