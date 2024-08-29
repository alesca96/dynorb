#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define GNCLIB_IMPLEMENTATION
#include "..\include\gnc101lib.h" // Include your custom RK header file

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

void SimpHarmOsc(const double t, const double yy[], void *params, double ff[])
{
    // // Parameters:
    // SimpHarmOscParams *p = (SimpHarmOscParams *)params;
    // double F0 = p->F0;
    // double m = p->m;
    // double om_n = p->om_n;
    // double zeta = p->zeta;
    // double om = p->om;

    SimpHarmOscParams p = {0};
    double F0 = 1.0;
    double m = 1.0;
    double om_n = 1.0;
    double om = 0.4 * om_n;
    double zeta = 0.03;
    // Simple Harmonic Oscillator: ff = df/dt
    ff[0] = yy[1];
    ff[1] = (F0 / m) * sin(om * t) - 2 * zeta * om_n * yy[1] - (om_n * om_n) * yy[0];
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

    return 0; // Success
}

int main(void)
{
    // Step 0: Define Parameters of ODE system:
    SimpHarmOscParams p = {0};
    p.F0 = 1.0;
    p.m = 1.0;
    p.om_n = 1.0;
    p.om = 0.4 * p.om_n;
    p.zeta = 0.03;

    // Step 1: Initial conditions and Integration Range:
    double t0 = 0.0;              // Initial time
    double t1 = 110.0;            // Final time
    double x0 = 0.0;              // Initial Position
    double x_dot0 = 0.0;          // Initial Velocity
    double yy0[2] = {x0, x_dot0}; // Initial State
    double h = 0.5;               // Step size
    const int sys_size = 2;       // System size (2 ODEs)
    const int rk_order = 4;       // RK4 method
    int num_steps = (int)((t1 - t0) / h) + 1;

    // Allocate memory for time and solution arrays:
    double *tt = (double *)malloc(num_steps * sizeof(double));
    double *yyt = (double *)malloc(num_steps * sys_size * sizeof(double));

    // Step 2: Perform Integration using custom RK method:
    gnc_rk1to4(SimpHarmOsc, t0, t1, h, rk_order, yy0, sys_size, tt, yyt);

    // Step 3: Open file to store results
    FILE *outfile = fopen("./data/ex_01_18b.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Step 4: Loop over time steps and calculate analytical solution
    for (int i = 0; i < num_steps; i++)
    {
        double t = tt[i];
        double yy[2] = {yyt[i * sys_size], yyt[i * sys_size + 1]};

        // Analytical Solution
        double x_a = 0.0;
        int status_analytical = SimpHarmOscAnalyticalSolution(t, yy0, &x_a, &p);
        if (status_analytical != 0)
        {
            printf("Error: Analytical solution calculation failed at t = %f\n", t);
            break;
        }

        fprintf(outfile, "%f %f %f %f\n", t, yy[0], yy[1], x_a);
    }

    // Step 5: Free Memory and Close Data File:
    free(tt);
    free(yyt);
    fclose(outfile);

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command =
        "set terminal qt\n"
        "set title 'Example 18 Chapter 01: Simple Harmonic Oscillator using Custom RK4'\n"
        "set xlabel 'Time t [s]'\n"
        "set ylabel 'x(t) [m], v(t) [m/s], x_a(t) [m]'\n"
        "plot './data/ex_01_18b.txt' using 1:2 with lines title 'x(t)', "
        "'./data/ex_01_18b.txt' using 1:3 with lines title 'v(t)', "
        "'./data/ex_01_18b.txt' using 1:4 with points pt 7 ps 1 lc rgb 'black' title 'x_a(t)'\n";

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    return 0;
}
