#define DYNORB_IMPLEMENTATION
#define USE_FLOAT // USE_real
#define USE_CBLAS
#include "..\include\dynorb.h" // Include your custom RK header file

/* Function to print test information: */
void printTestInfo(const char *description, int test_idx)
{
    printf("\n--------------------------------------------\n");
    printf("\033[1m\033[31mTest %d: %s\033[0m\n", test_idx, description);
}

/* Unit Test for all functions: */

int main(void)
{
    printf("\n");
    printf("\033[1m\033[34m==============================================\033[0m\n");
    printf("\033[1m\033[34mTEST SCRIPT begins:\033[0m\n");

    // Test Index:
    int test_number = 0;

    /* TEST: */
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

    /* TEST: */
    {
        const char *test_description = "Test Printing Matrix";
        printTestInfo(test_description, ++test_number);

        // Define a 2x3 matrix as a 1D array
        real matrix[12] = {1.0, 2.0, 3.0, 4.0, // transpose representation
                           5.0, 6.0, 7.0, 8.0,
                           9.0, 10.0, 11.0, 12.0};

        // Use the macro to access elements
        _dynorb_rmprint(matrix, 4, 3);
    }

    /* TEST: */
    {
        const char *test_description = "Test Copying Matrix";
        printTestInfo(test_description, ++test_number);

        int src_rows = 4, src_cols = 4;
        real src_M[16] = {
            1, 5, 9, 13,  // Column 1
            2, 6, 10, 14, // Column 2
            3, 7, 11, 15, // Column 3
            4, 8, 12, 16  // Column 4
        };

        int dst_rows = 4, dst_cols = 4;
        real dst_M[16] = {0}; // Destination initialized to all zeros

        // tCopy matrix starting at position (0,0) in the destination matrix
        _dynorb_rmcopy(src_M, src_rows, src_cols,
                       0, 0, 4, 4, // Submatrix starts at (1,1) and is 2x2
                       dst_M, dst_rows, dst_cols,
                       0, 0); // Copy it to (0,0) in the destination matrix

        // Use the macro to access elements
        printf("Source: \n");
        _dynorb_rmprint(src_M, src_rows, src_cols);
        printf("Destination: \n");
        _dynorb_rmprint(dst_M, dst_rows, dst_cols);

        printf("Element by element: \n");
        for (int i = 0; i < src_rows; ++i)
        {
            for (int j = 0; j < src_cols; ++j)
            {
                if (fabs(EL(src_M, src_rows, src_cols, i, j) - EL(dst_M, src_rows, src_cols, i, j)) < 1.e-12)
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

    /* TEST: */
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
        _dynorb_rmprint(C, 5, 1);
        printf("\n");
        printf("Resulting ROW vector C = alpha* A + B:\n");
        _dynorb_rmprint(C, 1, 5);
    }

    /* TEST: */
    {
        const char *test_description = "Test Cblas-like Matrix Vector Product - No Transposition";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (2x3) in column-major order
        real A[6] = {
            1.0, 2.0, // First column
            3.0, 4.0, // Second column
            5.0, 6.0  // Third column
        };

        // Define vectors X and Y
        real X[3] = {1.0, 1.0, 1.0}; // Vector X (size 3)
        real Y[2] = {0.0, 0.0};      // Vector Y (size 2)

        // Define scalars alpha and beta
        real alpha = 1.0;
        real beta = 0.0;

        // Perform matrix-vector multiplication Y = alpha * A * X + beta * Y
        _dynorb_rgemv(false, 2, 3, alpha, A, X, beta, Y);

        // Print the result
        printf("Y = [%f, %f]\n", Y[0], Y[1]);
    }

    /* TEST: */
    {
        const char *test_description = "Test Cblas-like Matrix Vector Product - Transposition";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (3x2) in column-major order
        real A[6] = {
            1.5, 4.7, // First column
            2.9, 5.6, // Second column
            3.3, 6.2  // Third column
        };

        // Define vectors X and Y
        real X[3] = {1.5, -6.6, 1.1}; // Vector X (size 3)
        real Y[2] = {0.0, 0.0};       // Vector Y (size 2)

        // Define scalars alpha and beta
        real alpha = 1.0;
        real beta = 0.0;

        // Perform matrix-vector multiplication Y = alpha * A^T * X + beta * Y
        _dynorb_rgemv(true, 3, 2, alpha, A, X, beta, Y);

        // Print the result
        printf("Y = [%f, %f]\n", Y[0], Y[1]);
    }

    /* TEST: Test Complted: */

    printf("\n--------------------------------------------\n");
    printf("\033[1m\033[32mTEST SCRIPT ends:\n");
    printf("ALL TEST COMPLETED without crashes.\n");
    printf("==============================================\033[0m\n");
    printf("\n");

    return 0;
}
