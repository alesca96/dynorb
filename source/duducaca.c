#include <stdio.h>

int main(void)
{
    int n = 10;
    printf("\nn = %d\n", n);
    int y[n];

    for (int i = 0; i < n; ++i)
    {
        y[i] = i;
        printf("y[%d] = %d\n", i, y[i]);
    }

    n = 12;
    printf("\nn = %d\n", n);

    int x[n];

    for (int i = 0; i < n; ++i)
    {
        x[i] = i;
        printf("x[%d] = %d\n", i, x[i]);
    }
}