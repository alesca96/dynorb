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
    // Parameters:
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

int SimpHarmOscAnalyticalSolution(double t, const double yy0[], double *x_analytical, void *params)
{

    // Initial Conditions:
    double x0 = yy0[0];
    double x_dot0 = yy0[1];
    // Parameters:
    SimpHarmOscParams *p = (SimpHarmOscParams *)params;
    double F0 = p->F0;
    double m = p->m;
    double om_n = p->om_n;
    double zeta = p->zeta;
    double om = p->om;
    // Intermediate Variables:
    double zeta2 = zeta * zeta;
    double om2 = om * om;
    double om_n2 = om_n * om_n;
    double omom_n = om * om_n;
    double om_d = om_n * sqrt(1 - zeta2);
    double _2omom_nzeta = (2 * omom_n * zeta);
    double F0m = (F0 / m);
    // Coefficients:
    double den = ((om_n2 - om2) * (om_n2 - om2)) + (_2omom_nzeta * _2omom_nzeta);
    double A = (zeta * (om_n / om_d) * x0) + (x_dot0 / om_d) + (((om2 + ((2 * zeta2 - 1) * om_n2)) / (den)) * (om / om_d) * F0m);
    double B = (x0) + ((_2omom_nzeta / den) * F0m);
    // Position x(t):
    *x_analytical = (exp(-zeta * om_n * t) * (A * sin(om_d * t) + B * cos(om_d * t))) + ((F0m / den) * (((om_n2 - om2) * sin(om * t)) - (_2omom_nzeta * cos(om * t))));

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
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, h, abstol, reltol);

    // Step 4: Open file to store results
    FILE *outfile = fopen("./data/ex_01_18-a.txt", "w");
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

        // Analytical Solution
        double x_a = 0.0;
        double yy0[2] = {x0, x_dot0};
        int status_analytical = SimpHarmOscAnalyticalSolution(t, yy0, &x_a, &p);
        if (status_analytical != GSL_SUCCESS)
        {
            printf("Error: %s\n", gsl_strerror(status));
            break;
        }

        fprintf(outfile, "%.10f %.10f %.10f %.10f\n", t, yy[0], yy[1], x_a);
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
        "set ylabel 'x(t) [m], v(t) [m/s], x_a(t) [m]'\n"
        "plot './data/ex_01_18-a.txt' using 1:2 with lines title 'x(t)', "
        "'./data/ex_01_18-a.txt' using 1:3 with lines title 'v(t)', "
        "'./data/ex_01_18-a.txt' using 1:4 with points pt 7 ps 1 lc rgb 'black' title 'x_a(t)'\n";

    // "plot './data/ex_01_18-a.txt' using 1:2 with points pt 1 ps 3 lc rgb 'red' title 'x(t)', "
    // "'./data/ex_01_18-a.txt' using 1:4 with points pt 2 ps 3 lc rgb 'blue' title 'x_a(t)'\n";

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    return 0;
}
