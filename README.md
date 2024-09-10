# dynorb
Implementation of Orbital, Attitude and GNC functions in C

![Logo](https://github.com/alesca96/dynorb/blob/main/logo_v00.png)

## ASSUMPTIONS
 * Code uses **COLUMN MAJOR** convention (see CBLAS).
 * Strides other than 1 are **not supported** in this implementation.
 * Type `real` can only be `float` (single precision) or `double` (double precision).

## TODO: Short/Medium Term
 * Implement Curtis Algorithms using wrappers of CBLAS and LAPACKE.
 * Reduce (completely eliminate when possible) heap usage.
 * Fulfill as much as possible NASA POWER OF TEN.

## TODO: Long Term
 * Implement Custom subroutines to substitute CBLAS and LAPACKE.
 * Completely eliminate heap usage (if impossible, remove the function from the library).
 * Enforce NASA POWER OF TEN.

## NASA POWER OF TEN RULES
(http://web.eecs.umich.edu/~imarkov/10rules.pdf):
 * Avoid complex flow constructs, such as goto and recursion.
 * All loops must have fixed bounds. This prevents runaway code.
 * Avoid heap memory allocation.
 * Restrict functions to a single printed page.
 * Use a minimum of two runtime assertions per function.
 * Restrict the scope of data to the smallest possible.
 * Check the return value of all non-void functions, or cast to void to indicate the return value is useless.
 * Use the preprocessor sparingly.
 * Limit pointer use to a single dereference, and do not use function pointers.
 * Compile with all possible warnings active; all warnings should then be addressed before release of the software.
