#define DYNORB_IMPLEMENTATION
#define USE_DOUBLE // #define USE_FLOAT //#define USE_DOUBLE
#define USE_CBLAS
#include "..\include\dynorb.h" // Include your custom RK header file

/* Function to print test information: */
void printTestInfo(const char *description, int test_idx)
{
    printf("\n--------------------------------------------\n");
    printf("\033[1m\033[31mTest %d: %s\033[0m\n\n", test_idx, description);
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
        const char *test_description = "Test MACRO: _dynorb_EL";
        printTestInfo(test_description, ++test_number);

        // Define a 3x2 matrix as a 1D array
        real A[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        // Use the macro to access elements
        printf("Element at (0, 0): %.3f\n", _dynorb_EL(A, 3, 2, 0, 0)); // 1
        printf("Element at (1, 1): %.3f\n", _dynorb_EL(A, 3, 2, 1, 0)); // 2
        printf("Element at (0, 1): %.3f\n", _dynorb_EL(A, 3, 2, 0, 1)); // 4
        printf("Element at (1, 2): %.3f\n", _dynorb_EL(A, 3, 2, 2, 1)); // 6
    }

    /* TEST: */
    {
        const char *test_description = "Test Simple Utility Functions";
        printTestInfo(test_description, ++test_number);

        int a = 5;
        int b = 6;
        real _a = 5.0;
        real _b = 6.0;
        printf("Values: (int) a = %d, (int) b = %d\n", a, b);
        printf("Values: (real) _a = %f, (real) _b = %f\n", _a, _b);
        printf("Maximum Integer: %d\n", _dynorb_max(a, b));
        printf("Minimum Integer: %d\n", _dynorb_min(a, b));
        printf("Maximum Real: %f\n", _dynorb_rmax(_a, _b));
        printf("Minimum Real: %f\n", _dynorb_rmin(_a, _b));
    }

    /* TEST: */
    {
        const char *test_description = "Test Printing Matrix";
        printTestInfo(test_description, ++test_number);

        // A [4x3]
        real A[4 * 3] = {1.0, 2.0, 3.0, 4.0,     // 1st column
                         5.0, 6.0, 7.0, 8.0,     // 2nd column
                         9.0, 10.0, 11.0, 12.0}; // 3rd column

        // B [5x1]
        real B[5 * 1] = {1.1, 1.2, 1.3, 1.4, 1.5};
        // C [1x5]
        real C[1 * 5] = {1.1, 1.2, 1.3, 1.4, 1.5};
        // Use the macro to access elements
        printf("Matrix A: [m=4, n=3]\n");
        _dynorb_rmprint(A, 4, 3);
        printf("Matrix B: [m=5, n=1]\n");
        _dynorb_rmprint(B, 5, 1);
        printf("Matrix C: [m=1, n=5]\n");
        _dynorb_rmprint(C, 1, 5);
    }

    /* TEST: */
    {
        const char *test_description = "Test Matrix Copying (uses Vector Copy)";
        printTestInfo(test_description, ++test_number);

        // Source:
        int src_rows = 4, src_cols = 4;
        real src_M[16] = {
            1, 5, 9, 13,  // Column 1
            2, 6, 10, 14, // Column 2
            3, 7, 11, 15, // Column 3
            4, 8, 12, 16  // Column 4
        };
        // Destination-full:
        int dst_rows1 = 4, dst_cols1 = 4;
        real dst_M1[16] = {0}; // Destination initialized to all zeros

        // Destination-bigger:
        int dst_rows2 = 5;
        int dst_cols2 = 6;
        real dst_M2[30] = {0};

        // Copy full source matrix starting at position (0,0) in the destination matrix
        _dynorb_rmcopy(src_M, src_rows, src_cols,    // Sourse shape
                       0, 0, 4, 4,                   // Copied (subportion) shape
                       dst_M1, dst_rows1, dst_cols1, // Destination shape
                       0, 0);                        // Starting point on Destination

        // Copy full source matrix starting at position (1,1) in the destination matrix
        _dynorb_rmcopy(src_M, src_rows, src_cols,    // Sourse shape
                       1, 1,                         // Starting point on Source
                       3, 3,                         // Copied (subportion) shape
                       dst_M2, dst_rows2, dst_cols2, // Destination shape
                       2, 2);                        // Starting point on Destination

        // Use the macro to access elements
        printf("Source Matrix:\n");
        _dynorb_rmprint(src_M, src_rows, src_cols);
        printf("\n");
        printf("Destination 1: Exact Copy of Source of same shape:\n");
        _dynorb_rmprint(dst_M1, dst_rows1, dst_cols1);
        printf("\n");
        printf("Destination 1: element by element analysis:\n");
        for (int i = 0; i < src_rows; ++i)
        {
            for (int j = 0; j < src_cols; ++j)
            {
                if (fabs(_dynorb_EL(src_M, src_rows, src_cols, i, j) - _dynorb_EL(dst_M1, src_rows, src_cols, i, j)) < 1.e-12)
                {
                    printf("Elements at position (%d, %d) are equal.\n", i, j);
                }
                else
                {
                    printf("Elements at position (%d, %d) are NOT equal.\n", i, j);
                }
            }
        }
        printf("\n");
        printf("Destination 2: Partial Copy of Source into bigger Destination:\n");
        _dynorb_rmprint(dst_M2, dst_rows2, dst_cols2);
    }

    /* TEST: */
    {
        const char *test_description = "Test Dot Product and Norm-2: ";
        printTestInfo(test_description, ++test_number);

        // Vectors:
        const int n = 3;
        real xx[3] = {1.0, 2.0, 3.0};
        real yy[3] = {4.0, 5.0, 6.0};

        // Results:
        real dot = _dynorb_rdot(n, xx, yy);
        real norm2_xx = _dynorb_rnrm2(n, xx);
        real norm2_yy = _dynorb_rnrm2(n, yy);

        // Print:
        printf("Vector xx (row representation):\n");
        _dynorb_rmprint(xx, 1, n);
        printf("Vector yy (row representation):\n");
        _dynorb_rmprint(yy, 1, n);
        printf("\n");
        printf("Dot Product: %f\n", dot);
        printf("Norm-2 (xx): %f\n", norm2_xx);
        printf("Norm-2 (xx): %f\n", norm2_yy);
    }

    /* TEST: */
    {
        const char *test_description = "Test Scaling of Vector (array actually): ";
        printTestInfo(test_description, ++test_number);

        // Scalar:
        real a = 3.0;

        // Dimensions:
        const int m = 3;
        const int n = 4;

        // Vector:
        real xx[4] = {1.0, 2.0, 1.5, 2.5};

        // Matrix:
        real A[3 * 4] = {
            1.0, 5.0, 9.0, 13.0,  // Column 1
            2.0, 6.0, 10.0, 14.0, // Column 2
            3.0, 7.0, 11.0, 15.0  // Column 3
        };

        // Print Initial Vector:
        printf("Initial Vector:\n");
        _dynorb_rmprint(xx, 1, n);
        printf("\n");
        printf("Initial Matrix:\n");
        _dynorb_rmprint(A, m, n);
        printf("\n");
        // Scale Vector and Matrix:
        _dynorb_rscal(n, a, xx);
        _dynorb_rscal(m * n, a, A);
        printf("Scaled Vector by a = %f:\n", a);
        _dynorb_rmprint(xx, 1, n);
        printf("\n");
        printf("Scaled Matrix by a = %f:\n", a);
        _dynorb_rmprint(A, m, n);
        printf("\n");
    }

    /* TEST: */
    {
        const char *test_description = "Test Vector Sum";
        printTestInfo(test_description, ++test_number);

        // Vectors:
        int n = 5;
        real xx[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
        real yy[5] = {5.0, 4.0, 3.0, 2.0, 1.0};
        real zz1[5] = {0};
        real zz2[5] = {0};
        memcpy(zz1, yy, 5 * sizeof(real));
        memcpy(zz2, yy, 5 * sizeof(real));

        // Compute zz1 = 1.0*xx + zz1 with (zz1=yy)
        _dynorb_raxpy(n, 1.0, xx, zz1);
        // Compute zz2 = 3.5*xx + zz2 with (zz2=yy)
        _dynorb_raxpy(n, 3.5, xx, zz2);

        // Print the result
        printf("Resulting vector zz1 = 1.0*xx + zz1, printed as column:\n");
        _dynorb_rmprint(zz1, 5, 1);
        printf("\n");
        printf("Resulting vector zz2 = 3.5*xx + zz2, printed as row:\n");
        _dynorb_rmprint(zz2, 1, 5);
        printf("\n");
    }

    /* TEST: */
    {
        const char *test_description = "Test Matrix Vector Product - No Transposition";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (2x3) in column-major order
        real A[6] = {1.5, 4.7,
                     2.9, 5.6,
                     3.3, 6.2};

        // Define vectors X and Y
        real xx[3] = {1.5, -6.6, 1.1}; // Vector X (size 3)
        real yy[2] = {0.0, 0.0};       // Vector Y (size 2)

        // Define scalars alpha and beta
        real alpha = 1.0;
        real beta = 0.0;

        // Perform matrix-vector multiplication Y = alpha * A * X + beta * Y
        _dynorb_rgemv(false, 2, 3, alpha, A, xx, beta, yy);

        // Print the result
        printf("Y = [%f, %f]\n", yy[0], yy[1]);
    }

    /* TEST: */
    {
        const char *test_description = "Test Matrix Vector Product - Transposition";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (2x3) in column-major order
        real A[6] = {1.5, 4.7,
                     2.9, 5.6,
                     3.3, 6.2};

        // Define vectors X and Y
        real xx[2] = {2.3, 6.7};      // Vector X (size 2)
        real yy[3] = {0.0, 0.0, 0.0}; // Vector Y (size 3)

        // Define scalars alpha and beta
        real alpha = 1.0;
        real beta = 0.0;

        // Perform matrix-vector multiplication Y = alpha * A * X + beta * Y
        _dynorb_rgemv(true, 2, 3, alpha, A, xx, beta, yy);

        // Print the result
        printf("Y = [%f, %f, %f]\n", yy[0], yy[1], yy[2]);
    }

    /* TEST: */
    {
        const char *test_description = "Test Matrix Vector Product - Transposition - alpha 3.0 beta 1.0";
        printTestInfo(test_description, ++test_number);

        // Define matrix A (2x3) in column-major order
        real A[6] = {1.5, 4.7,
                     2.9, 5.6,
                     3.3, 6.2};

        // Define vectors X and Y
        real xx[2] = {2.3, 6.7};      // Vector X (size 2)
        real yy[3] = {0.0, 0.0, 0.0}; // Vector Y (size 3)

        // Define scalars alpha and beta
        real alpha = 3.0;
        real beta = 1.0;

        // Perform matrix-vector multiplication Y = alpha * A * X + beta * Y
        _dynorb_rgemv(true, 2, 3, alpha, A, xx, beta, yy);

        // Print the result
        printf("Y = [%f, %f, %f]\n", yy[0], yy[1], yy[2]);
    }

    /* TEST: */
    {
        const char *test_description = "Test Static Memory Allocation and Fill Vector";
        printTestInfo(test_description, ++test_number);

        // Initialize:
        _dynorb_solverConf solv_conf = {
            .h = 1.0,
            .n_steps = (int)(((10.0 - 0.0) + solv_conf.h - 1.0) / solv_conf.h),
        };

        // Use the calculated n_steps to declare the array
        real xx[solv_conf.n_steps];
        // Fill vector with zeros:
        _dynorb_rvfill(xx, (solv_conf.n_steps), 0.0);
        // Print
        printf("xx =\n");
        _dynorb_rmprint(xx, 1, (solv_conf.n_steps));

        // Use the calculated n_steps to declare the array
        const int n_steps = solv_conf.n_steps;
        real yy[n_steps];
        // Fill other values:
        _dynorb_rvfill(yy, n_steps, 42.0);
        // Print:
        printf("\nyy =\n");
        _dynorb_rmprint(yy, (solv_conf.n_steps), 1);
    }

    /* TEST:  Complted: */
    printf("\n--------------------------------------------\n");
    printf("\033[1m\033[32mTEST SCRIPT ends:\n");
    printf("ALL TEST COMPLETED without crashes.\n");
    printf("==============================================\033[0m\n");
    printf("\n");
    return 0;
}
