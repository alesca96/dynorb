#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

/* Define the system of ODEs: example 1.18 (Chapter 1, pag.45)
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.) */

int SimpHarmOsc(double x, const double y[], double f[], void *params)
{
    // Simple Harmonic Oscillator:
    f[0] = y[1];
    f[1] = -y[0];

    return GSL_SUCCESS;
}

int main(void)
{

    /* GSL: Use GSL to Integrate the ODE system: */

    // Step 1: Define the ODE sys:
    gsl_odeiv2_system sys = {SimpHarmOsc, NULL, 2, NULL};

    // Step 2: Choose a Solver (RK4) and setup stepping function:
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1.e-6, 1.e-6, 0.0);

    // Step 3: Initial conditions and Integration Range:
    double x = 0.0;           // Initial x
    double x1 = 10.0;         // End x
    double y[2] = {1.0, 0.0}; // Initial y1 and y2
    double h = 1e-6;          // Step size

    // Step 4: Open file to store results
    FILE *outfile = fopen("./data/ex_01_18.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Step 5: Perform Integration:
    for (int i = 0; i <= 100; ++i)
    {
        double xi = i * (x1 / 100.0);
        int status = gsl_odeiv2_driver_apply(driver, &x, xi, y);
        if (status != GSL_SUCCESS)
        {
            printf("Error: %s\n", gsl_strerror(status));
            break;
        }
        // printf("x = %.5f, y1 = %.5f, y2 = %.5f\n", x, y[0], y[1]); // Print results
        fprintf(outfile, "%f %f %f\n", x, y[0], y[1]);
    }

    // Step 6: Free Memory and Close Data File:
    gsl_odeiv2_driver_free(driver);
    fclose(outfile);

    /* GNUPLOT: Use Gnuplot to plot the data */
    const char *plot_command =
        "set terminal qt\n"
        "set title 'Example 18 Chapter 01: Simple Harmonic Oscillator using GSL'\n"
        "set xlabel 'Time, t'\n"
        "set ylabel 'x(t), v(t)'\n"
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
