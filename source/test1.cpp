#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

// Define the system of ODEs: dy/dx = -y
int func(double x, const double y[], double f[], void *params)
{
    f[0] = -y[0];
    return GSL_SUCCESS;
}

int main()
{
    printf("===============================\n");
    printf("Program Start\n");
    printf("-------------------------------\n");

    // Initial conditions
    double y[1] = {1.0}; // Initial value of y
    double x = 0.0;      // Initial x
    double x_end = 5.0;  // End x
    double h = 1e-6;     // Step size

    // GSL ODE system setup
    gsl_odeiv2_system sys = {func, NULL, 1, NULL};

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-6, 1e-6);

    // Open file to store results
    FILE *outfile = fopen("./data/data.txt", "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Integrate the ODE
    for (int i = 0; i <= 100; i++)
    {
        double xi = i * x_end / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &x, xi, y);

        if (status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error: %s\n", gsl_strerror(status));
            break;
        }

        // Store the results for plotting
        fprintf(outfile, "%f %f\n", x, y[0]);
        // printf("x = %f, y = %f\n", x, y[0]);
    }

    gsl_odeiv2_driver_free(d);
    fclose(outfile);

    // Use Gnuplot to plot the data
    const char *plot_command =
        "set terminal qt\n"
        "set title 'Solution of dy/dx = -y using GSL'\n"
        "set xlabel 'x'\n"
        "set ylabel 'y'\n"
        "plot './data/data.txt' with lines title 'y(x)'\n";

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    printf("-------------------------------\n");
    printf("Program End\n");
    printf("===============================\n");

    return 0;
}
