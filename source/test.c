#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h" // Include your custom RK header file

// Function to print test information
void printTestInfo(const char *description, int test_idx)
{
    printf("\n--------------------------------------------\n");
    printf("Test %d: %s\n", test_idx, description);
}

int main(void)
{
    printf("\n");
    printf("==============================================\n");
    printf("TEST SCRIPT begins:\n");

    // Test Index:
    int test_number = 0;

    /* TEST 0: */
    {
        const char *test_description = "Test Cblas-like flags";
        printTestInfo(test_description, test_number);
        printf("CblasColMajor = %d | CblasNoTrans = %d\n", CblasColMajor, CblasNoTrans);
    }

    /* TEST 1: */
    {
        const char *test_description = "Test macro EL";
        printTestInfo(test_description, ++test_number);

        // Define a 3x2 matrix as a 1D array
        real matrix[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        // Use the macro to access elements
        printf("Element at (0, 0): %.3f\n", EL(matrix, 3, 2, 0, 0)); // 1
        printf("Element at (1, 1): %.3f\n", EL(matrix, 3, 2, 1, 0)); // 2
        printf("Element at (0, 1): %.3f\n", EL(matrix, 3, 2, 0, 1)); // 4
        printf("Element at (1, 2): %.3f\n", EL(matrix, 3, 2, 2, 1)); // 6
    }

    /* TEST 2: */
    {
        const char *test_description = "Test Printing Matrix";
        printTestInfo(test_description, ++test_number);

        // Define a 2x3 matrix as a 1D array
        real matrix[12] = {1.0, 2.0, 3.0, 4.0, // transpose representation
                           5.0, 6.0, 7.0, 8.0,
                           9.0, 10.0, 11.0, 12.0};

        // Use the macro to access elements
        matprint(matrix, 4, 3);

        matprint(matrix, 3, 4);
    }

    /* TEST 3: */
    {
        const char *test_description = "Test Copying Matrix";
        printTestInfo(test_description, ++test_number);

        // Define a 2x3 matrix as a 1D array
        real src_M[12] = {1.0, 2.0, 3.0, 4.0, // transpose representation
                          5.0, 6.0, 7.0, 8.0,
                          9.0, 10.0, 11.0, 12.0};
        real dst_M[12] = {0.0};

        // Use the macro to access elements
        int nrows = 4;
        int ncols = 3;
        printf("Source: \n");
        matprint(src_M, nrows, ncols);
        printf("Destination: \n");
        matprint(dst_M, nrows, ncols);
        printf("Element by element: \n");
        for (int i = 0; i < nrows; ++i)
        {
            for (int j = 0; j < ncols; ++j)
            {
                if (fabs(EL(src_M, nrows, ncols, i, j) - EL(src_M, nrows, ncols, i, j)) < 1.e-12)
                {
                    printf("Elements at position (%d, %d) are equal.\n", i, j);
                }
                else
                {
                    printf("Elements at position (%d, %d) are NOT equal.\n", i, j);
                }
            }
        }
    }

    /* TEST 4: */
    {
        const char *test_description = "Test Cblas-like Vector Sum";
        printTestInfo(test_description, ++test_number);

        // Define the vectors and their size
        int n = 5;
        real A[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
        real B[5] = {5.0, 4.0, 3.0, 2.0, 1.0};
        real C[5];
        memcpy(C, B, 5 * sizeof(real));

        // Use cblas_saxpy to perform C = A + C (where C = B initially)
        _dynorb_raxpy(n, 1.0, A, C);

        // Print the result
        printf("Resulting COLUMN vector C = alpha* A + B:\n");
        matprint(C, 5, 1);
        printf("\n");
        printf("Resulting ROW vector C = alpha* A + B:\n");
        matprint(C, 1, 5);
    }

    /* TEST 6: */
    {
        const char *test_description = "Test Cblas-like Matrix Vector Product - No Transposition";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (2x3) in column-major order
        double A[6] = {
            1.0, 2.0, // First column
            3.0, 4.0, // Second column
            5.0, 6.0  // Third column
        };

        // Define vectors X and Y
        double X[3] = {1.0, 1.0, 1.0}; // Vector X (size 3)
        double Y[2] = {0.0, 0.0};      // Vector Y (size 2)

        // Define scalars alpha and beta
        double alpha = 1.0;
        double beta = 0.0;

        // Perform matrix-vector multiplication Y = alpha * A * X + beta * Y
        _dynorb_rgemv(false, 2, 3, alpha, A, X, beta, Y);

        // Print the result
        printf("Y = [%f, %f]\n", Y[0], Y[1]);
    }

    /* TEST 7: */
    {
        const char *test_description = "Test Cblas-like Matrix Vector Product - Transposition";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (3x2) in column-major order
        double A[6] = {
            1.5, 4.7, // First column
            2.9, 5.6, // Second column
            3.3, 6.2  // Third column
        };

        // Define vectors X and Y
        double X[3] = {1.5, -6.6, 1.1}; // Vector X (size 3)
        double Y[2] = {0.0, 0.0};       // Vector Y (size 2)

        // Define scalars alpha and beta
        double alpha = 1.0;
        double beta = 0.0;

        // Perform matrix-vector multiplication Y = alpha * A^T * X + beta * Y
        _dynorb_rgemv(true, 3, 2, alpha, A, X, beta, Y);

        // Print the result
        printf("Y = [%f, %f]\n", Y[0], Y[1]);
    }

    /* TEST N: Test Complted: */

    printf("\n--------------------------------------------\n");
    printf("TEST SCRIPT ends:\n");
    printf("ALL TEST COMPLETED without crashes \n");
    printf("==============================================\n");
    printf("\n");

    return 0;
}
