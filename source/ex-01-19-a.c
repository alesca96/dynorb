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

void SimpHarmOsc(const real in_t, const real *in_yy, const void *in_params, real *out_ff)
{
    // Parameters:
    SimpHarmOscParams *p = (SimpHarmOscParams *)in_params;
    real F0 = p->F0;
    real m = p->m;
    real om_n = p->om_n;
    real zeta = p->zeta;
    real om = p->om;

    // Simple Harmonic Oscillator: out_ff = dyy/dt
    out_ff[0] = in_yy[1];
    out_ff[1] = (F0 / m) * sin(om * in_t) - 2 * zeta * om_n * in_yy[1] - (om_n * om_n) * in_yy[0];
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
    // Step 0: Define Parameters of ODE system:
    SimpHarmOscParams p =
        {
            .F0 = 1.0,
            .m = 1.0,
            .om_n = 1.0,
            .om = 0.4 * p.om_n,
            .zeta = 0.03,
        };

    // Step 1: Initial conditions and Step Size:
    const real x0 = 0.0;              // Initial Position
    const real x_dot0 = 0.0;          // Initial Velocity
    const real yy0[2] = {x0, x_dot0}; // Initial State

    // Step 2: Setup odeSys and solverConf structure:
    _dynorb_odeSys SimpHarmOscSys =
        {
            .odeFunction = (_dynorb_odeFun *)SimpHarmOsc,
            .params = &p,
            .sys_size = 2,
            .t0 = 0.0,
            .t1 = 110.0,
            .yy0 = yy0,
            .tt = NULL,
            .YY_t = NULL,
        };

    _dynorb_solverConf solv_conf =
        {
            .h = 0.1,
            .n_steps = (int)((((SimpHarmOscSys.t1 - SimpHarmOscSys.t0)) / solv_conf.h) + 1.0),
        };

    // Memory Allocation For Solution:
    real tt[solv_conf.n_steps];
    real Y_tt[solv_conf.n_steps * SimpHarmOscSys.sys_size];
    // Add Pointers to structure:
    SimpHarmOscSys.tt = tt;
    SimpHarmOscSys.YY_t = Y_tt;

    // Print:
    printf("Time Step Size: <h = %f [s]>\n", solv_conf.h);
    printf("Number Steps: <solv_conf.n_steps = %d>\n", solv_conf.n_steps);

    // Step 3: Perform Integration using custom RK method:
    _dynorb_heun_(&SimpHarmOscSys, &solv_conf);

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

    // Closing the file:
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
