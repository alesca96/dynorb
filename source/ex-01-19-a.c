#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
// #define USE_CBLAS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>            // For memcpy
#include "..\include\dynorb.h" // Include your custom RK header file

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

void SimpHarmOsc(const double in_t, const double *in_yy, const void *in_params, double *out_ff)
{
    // Parameters:
    SimpHarmOscParams *p = (SimpHarmOscParams *)in_params;
    double F0 = p->F0;
    double m = p->m;
    double om_n = p->om_n;
    double zeta = p->zeta;
    double om = p->om;

    // Simple Harmonic Oscillator: out_ff = dyy/dt
    out_ff[0] = in_yy[1];
    out_ff[1] = (F0 / m) * sin(om * in_t) - 2 * zeta * om_n * in_yy[1] - (om_n * om_n) * in_yy[0];
}

void SimpHarmOscAnalyticalSolution(double t, const double *yy0, double *x_analytical, const void *params)
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
}

int main(void)
{
    // Step 0: Define Parameters of ODE system:
    SimpHarmOscParams p =
        {
            .F0 = 1.0,
            .m = 1.0,
            .om_n = 1.0,
            .om = 0.4 * p.om_n,
            .zeta = 0.03

        };

    // Step 1: Initial conditions:
    const int sys_size = 2;             // System size (2 ODEs)
    const double t0 = 0.0;              // Initial time
    const double t1 = 110.0;            // Final time
    const double x0 = 0.0;              // Initial Position
    const double x_dot0 = 0.0;          // Initial Velocity
    const double yy0[2] = {x0, x_dot0}; // Initial State

    // Step 2: Collect Data into _dynorb_odeSys structure:
    _dynorb_odeSys SimpHarmOscSys = {
        .odeFunction = (const _dynorb_odeFun *)SimpHarmOsc,
        .params = &p,
        .sys_size = sys_size,
        .t0 = t0,
        .t1 = t1,
        .yy0 = yy0};

    // Step 3: Set up Integration:
    double h = 1.0; // 1.0; // Step size
    int num_steps = (ceil((t1 - t0) / h) + 1);
    printf("Time Step Size: <h = %f [s]>\n", h);
    printf("Number Steps: <num_steps = %d>\n", num_steps);

    // return 0;
    double *tt = (double *)malloc(num_steps * sizeof(double));             // Allocate memory for time array:
    double *yyt = (double *)malloc(num_steps * sys_size * sizeof(double)); // Allocate memory for solution array:

    // Step 5: Perform Integration using custom RK method:
    _dynorb_heun_(&SimpHarmOscSys, h, num_steps, tt, yyt);

    // Step 5: Open file to store results
    FILE *outfile = fopen("./data/ex_01_19a.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Step 6: Loop over time steps and calculate analytical solution
    for (int i = 0; i < num_steps - 1; i++) // -1 beacause last line might be spurious
    {
        // Analytical Solution
        double x_a = 0.0;
        SimpHarmOscAnalyticalSolution(tt[i], yy0, &x_a, &p);
        fprintf(outfile, "%.10f %.10f %.10f %.10f\n", tt[i], yyt[i * sys_size], yyt[i * sys_size + 1], x_a);
    }

    // Closing the file:
    fclose(outfile);

    // Step 7: Free Memory and Close Data File:
    free(yyt);
    free(tt);

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command =
        "set terminal qt enhanced\n"
        "set title 'Example 19 Chapter 01: Simple Harmonic Oscillator using Custom Huen'\n"
        "set xlabel 'Time t [s]'\n"
        "set ylabel 'x(t) [m], x_{a}(t) [m]'\n"
        "plot './data/ex_01_19a.txt' using 1:2 with lines lc rgb 'red' title 'x(t)', "
        "'./data/ex_01_19a.txt' using 1:4 with lines lc rgb 'black' title 'x_{a}(t)'\n";

    // "'./data/ex_01_19a.txt' using 1:3 with points pt 7 ps 0.5 lc rgb 'blue' title 'v(t)', "
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    return 0;
}
