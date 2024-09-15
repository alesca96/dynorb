#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE // USE_FLOAT
#define USE_CBLAS
#include "..\include\dynorb.h" // Include your custom RK header file

/* Define the system of ODEs: example 1.18 (Chapter 1, pag.45)
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.) */

typedef struct
{
    real F0;
    real m;
    real om_n;
    real zeta;
    real om;

} SimpHarmOscParams;

void SimpHarmOsc(const real t, const real *yy, const void *params, real *ff)
{
    // Parameters:
    SimpHarmOscParams *p = (SimpHarmOscParams *)params;
    real F0 = p->F0;
    real m = p->m;
    real om_n = p->om_n;
    real zeta = p->zeta;
    real om = p->om;

    // Simple Harmonic Oscillator: ff = dyy/dt
    ff[0] = yy[1];
    ff[1] = (F0 / m) * sin(om * t) - 2 * zeta * om_n * yy[1] - (om_n * om_n) * yy[0];
}

void SimpHarmOscAnalyticalSolution(real t, const real *yy0, real *x_analytical, const void *params)
{

    // Initial Conditions:
    real x0 = yy0[0];
    real x_dot0 = yy0[1];
    // Parameters:
    SimpHarmOscParams *p = (SimpHarmOscParams *)params;
    real F0 = p->F0;
    real m = p->m;
    real om_n = p->om_n;
    real zeta = p->zeta;
    real om = p->om;
    // Intermediate Variables:
    real zeta2 = zeta * zeta;
    real om2 = om * om;
    real om_n2 = om_n * om_n;
    real omom_n = om * om_n;
    real om_d = om_n * sqrt(1 - zeta2);
    real _2omom_nzeta = (2 * omom_n * zeta);
    real F0m = (F0 / m);
    // Coefficients:
    real den = ((om_n2 - om2) * (om_n2 - om2)) + (_2omom_nzeta * _2omom_nzeta);
    real A = (zeta * (om_n / om_d) * x0) + (x_dot0 / om_d) + (((om2 + ((2 * zeta2 - 1) * om_n2)) / (den)) * (om / om_d) * F0m);
    real B = (x0) + ((_2omom_nzeta / den) * F0m);
    // Position x(t):
    *x_analytical = (exp(-zeta * om_n * t) * (A * sin(om_d * t) + B * cos(om_d * t))) + ((F0m / den) * (((om_n2 - om2) * sin(om * t)) - (_2omom_nzeta * cos(om * t))));
}

int main(void)
{
    /* ==========================================================
     * DYNORB: Numerical Integration
     * ========================================================== */
    // Step 1: Set up of ode-System and Solver Configuration:

    // Declare Structures:
    _dynorb_odeSys SimpHarmOscSys;
    _dynorb_solverConf solv_conf;

    // Structure Fields:
    SimpHarmOscParams p =
        {
            .F0 = 1.0,
            .m = 1.0,
            .om_n = 1.0,
            .om = 0.4 * p.om_n,
            .zeta = 0.03,
        };
    real x0 = 0.0;              // Initial Position
    real x_dot0 = 0.0;          // Initial Velocity
    real yy0[2] = {x0, x_dot0}; // Initial State
    int sys_size = 2;           // Size of odeSys
    real t0 = 0.0;              // Initial Time
    real t1 = 110.0;            // Final Time
    real h = 1.0;               // Time step

    // Confugure (dynamically)
    _dynorb_configure_dynamic(&SimpHarmOscSys, &solv_conf,
                              &SimpHarmOsc, &p, yy0,
                              sys_size, t0, t1, h);

    // Step 3: Perform Integration using custom RK method:
    _dynorb_heun_(&SimpHarmOscSys, &solv_conf);

    /* ==========================================================
     * DATA LOG: Save data
     * ========================================================== */

    // Step 4: Open file to store results
    FILE *outfile = fopen("./data/ex_01_19-a.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Step 5: Loop over time steps and calculate analytical solution
    for (int i = 0; i < solv_conf.n_steps; i++) // -1 beacause last line might be spurious
    {
        // Analytical Solution
        real x_a = 0.0;
        SimpHarmOscAnalyticalSolution(SimpHarmOscSys.tt[i], yy0, &x_a, &p);
        fprintf(outfile, "%.10f %.10f %.10f %.10f\n", SimpHarmOscSys.tt[i], SimpHarmOscSys.YY_t[i * SimpHarmOscSys.sys_size], SimpHarmOscSys.YY_t[i * SimpHarmOscSys.sys_size + 1], x_a);
    }

    // Free Memory and Closing the file:
    _dynorb_free(&SimpHarmOscSys);
    fclose(outfile);

    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command =
        "set terminal qt enhanced\n"
        "set title 'Example 19 Chapter 01: Simple Harmonic Oscillator using Custom Huen'\n"
        "set xlabel 'Time t [s]'\n"
        "set ylabel 'x(t) [m], x_{a}(t) [m]'\n"
        "plot './data/ex_01_19-a.txt' using 1:2 with lines lc rgb 'red' title 'x(t)', "
        "'./data/ex_01_19-a.txt' using 1:4 with lines lc rgb 'black' title 'x_{a}(t)'\n";

    // "'./data/ex_01_19-a.txt' using 1:3 with points pt 7 ps 0.5 lc rgb 'blue' title 'v(t)', "
    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    return 0;
}
