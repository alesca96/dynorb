#include <stdio.h>

int main(void)
{
    /* ==========================================================
     * GNUPLOT: Use Gnuplot to plot the data
     * ========================================================== */

    const char *plot_command =
        "set terminal qt\n"
        "set title 'Example 18 Chapter 01: Compare GSL to Custom RK4'\n"
        "set xlabel 'Time t [s]'\n"
        "set ylabel 'x_{gsl}(t),x_{custom}(t), x_{analytical}(t)'\n"
        "plot './data/ex_01_18a.txt' using 1:2 with points pt 1 ps 1.5 lc rgb 'green' title 'x_{gsl}(t)', "
        "'./data/ex_01_18b.txt' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'x_{custom}(t)', "
        "'./data/ex_01_18b.txt' using 1:4 with lines lc rgb 'red' title 'x_{analytical}(t)'\n";

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot)
    {
        fprintf(gnuplot, "%s", plot_command);
        pclose(gnuplot);
    }

    return 0;
}