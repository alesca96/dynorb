#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
// #define USE_CBLAS
#include "..\include\dynorb.h" // Include your custom RK header file

int main(void)
{
    printf("Hello world | CblasColMajor = %d | CblasNoTrans = %d\n", CblasColMajor, CblasNoTrans);

    // Define the vectors and their size
    int n = 5;
    real A[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    real B[5] = {5.0, 4.0, 3.0, 2.0, 1.0};

    // Use cblas_saxpy to perform C = A + C (where C = B initially)
    _dynorb_raxpy(n, 1.0, A, 1, B, -1);

    // Print the result
    printf("Resulting vector C = alpha* A + B:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%f ", B[i]);
    }
    printf("\n");

    return 0;
}
