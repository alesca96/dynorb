#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h"
// Time integration
#include <sys/time.h>

/*  Example 2.16 (Chapter 2, pag.129): Locate the five Lagrange points for the
    Earth-Moon system using the Bisection method. ex-02-16-a.c
    REF: Curtis, H.D., 2020. Orbital mechanics for engineering students (3rd Edit.)*/

typedef struct
{
    real massRatio_2;
} EarthMoonParams;

real LagrangeEarthMoon(real eps, void *params)
{
    EarthMoonParams *p = (EarthMoonParams *)params;
    real f = ((1.0 - p->massRatio_2) * (eps + p->massRatio_2) / pow(fabs(eps + p->massRatio_2), 3.0)) + (p->massRatio_2 * (eps + p->massRatio_2 - 1) / pow(fabs(eps + p->massRatio_2 - 1), 3.0)) - eps;
    return f;
}

int main(void)
{
    // Input data :
    real m1 = 5.974e24;             // Earth mass
    real m2 = 7.348e22;             // Moon mass
    real r12 = 3.844e5;             // Earth-Moon distance
    real xxl[3] = {-1.1, 0.5, 1.0}; // Array of lower bounds
    real xxu[3] = {-0.9, 1.0, 1.5}; // Array of upper bounds

    // Ratio Pi2
    EarthMoonParams p;
    p.massRatio_2 = m2 / (m1 + m2);

    printf("\n=================================================\n");
    // Bisection method:
    real L[3];
    for (int i = 0; i < 3; ++i)
    {
        L[i] = _dynorb_bisect(LagrangeEarthMoon, &p, xxl[i], xxu[i], 1.e-6);
        printf("Episolon Lagrange Point L(xxl = %g, xxu= %g) = %g\n", xxl[i], xxu[i], L[i]);
    }

    // Print out the x-coordinates of L1, L2 and L3 relative to the center of mass.

    // Output to the command window:
    printf("\n Parameters: \n");
    printf("m1 = %g kg\n", m1);
    printf("m2 = %g kg\n", m2);
    printf("r12 = %g km\n", r12);
    printf("\n The 3 collinear Lagrange points:");
    printf("\n L3 : x = % 10g km (f(x3) = % g)", L[0] * r12, LagrangeEarthMoon(L[0], &p));
    printf("\n L1 : x = % 10g km (f(x2) = % g)", L[1] * r12, LagrangeEarthMoon(L[1], &p));
    printf("\n L2 : x = % 10g km (f(x1) = % g)", L[2] * r12, LagrangeEarthMoon(L[2], &p));
    printf("\n=================================================\n");
}