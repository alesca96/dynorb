
void _dynorb_rrkf45(_dynorb_odeSys *sys, _dynorb_solverConf *solver_configuration, const real tol)
{
    // Print WARNING:
    printf("\nWARNING: Arrays sys->tt and sys->YY_t must be dynamically allocated to use this function.\n");
    printf("         Please make sure to configure using _dynorb_configure_dynamic.\n");
    printf("         If this is the case already, ignore this warning.\n\n");

    // Open up structures [TODO: remove these copies]
    _dynorb_odeFun *odeFunction = sys->odeFunction; // Pointer to _dynorb_odeFun
    const void *odeParams = sys->odeParams;         // Pointer to odeParams
    const int sys_size = sys->sys_size;             // Size of System
    const real *yy0 = sys->yy0;                     // Pointer to Initial Conditions
    real t0 = sys->t0;                              // Initial time
    real t1 = sys->t1;                              // Final time
    real *tt = sys->tt;                             // Time steps of solution
    real *YY_t = sys->YY_t;                         // Solution Array
    real h = solver_configuration->h;               // Initial step
    int n_steps = solver_configuration->n_steps;

    // Coefficients:
    real a[6] = {0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0};
    real b[36] = {
        0.0, 1.0 / 4.0, 3.0 / 32.0, 1932.0 / 2197.0, 439.0 / 216.0, -8.0 / 27.0, // first column
        0.0, 0.0, 9.0 / 32.0, -7200.0 / 2197.0, -8.0, 2.0,                       // second column
        0.0, 0.0, 0.0, 7296.0 / 2197.0, 680.0 / 513.0, -3544.0 / 2565.0,         // third column
        0.0, 0.0, 0.0, 0.0, -845.0 / 4104.0, 1859.0 / 4104.0,                    // fourth column
        0.0, 0.0, 0.0, 0.0, 0.0, -11.0 / 40.0,                                   // fifth column
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0                                             // sixth column
    };
    real c_star[6] = {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0};
    real c[6] = {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};

    // Allocate Memory for Current State, Inner State:
    real yy[sys_size]; // High order solution (RK5)
    // real yy_star[sys_size]; // Low order solution (RK4)
    real yyi_[sys_size];
    real yy_inner[sys_size];
    real ee[sys_size]; // Truncation error vector
    real ff1_[sys_size];
    real ff2_[sys_size];
    real ff3_[sys_size];
    real ff4_[sys_size];
    real ff5_[sys_size];
    real ff6_[sys_size];

    // Allocate memory for scalars:
    real t;
    real ti;
    real t_inner;
    real hmin = 0.0;
    real trunc_err_max;
    real trunc_err_allowed;
    real delta;
    int step = 0;

    // Initialize:
    t = t0;
    // h = (t1 - t0) / 100.0;
    _dynorb_rvcopy(sys_size, yy0, yy); // Current state at t0 = initial state

    while (t < t1)
    {
        // printf("\n\nBegin time step t = %f | h = %f\n", t, h);
        hmin = 16.0 * (_dynorb_eps(t)); // Update hmin
        ti = t;
        _dynorb_rvcopy(sys_size, yy, yyi_);
        // Evaluate time derivatives at stage 1:
        t_inner = ti;
        _dynorb_rvcopy(sys_size, yyi_, yy_inner); // yy_inner = yyi_
        odeFunction(t_inner, yy_inner, odeParams, ff1_);

        // Evaluate time derivatives at stage 2:
        t_inner = ti + (a[1] * h);
        _dynorb_rvcopy(sys_size, yyi_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 1, 0)), ff1_, yy_inner); // yy_inner = yyi_ + b(2,1)*ff1_
        odeFunction(t_inner, yy_inner, odeParams, ff2_);

        // Evaluate time derivatives at stage 3:
        t_inner = ti + (a[2] * h);
        _dynorb_rvcopy(sys_size, yyi_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 2, 0)), ff1_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 2, 1)), ff2_, yy_inner); // yy_inner = yyi_ + b(3,1)*ff1_ + b(3,2)*ff2_
        odeFunction(t_inner, yy_inner, odeParams, ff3_);

        // Evaluate time derivatives at stage 4:
        t_inner = ti + (a[3] * h);
        _dynorb_rvcopy(sys_size, yyi_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 3, 0)), ff1_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 3, 1)), ff2_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 3, 2)), ff3_, yy_inner); // yy_inner = yyi_ + b(4,1)*ff1_ + b(4,2)*ff2_ + b(4,3)*ff3_
        odeFunction(t_inner, yy_inner, odeParams, ff4_);

        // Evaluate time derivatives at stage 5:
        t_inner = ti + (a[4] * h);
        _dynorb_rvcopy(sys_size, yyi_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 4, 0)), ff1_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 4, 1)), ff2_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 4, 2)), ff3_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 4, 3)), ff4_, yy_inner); // yy_inner = yyi_ + b(5,1)*ff1_ + b(5,2)*ff2_ + b(5,3)*ff3_ + b(5,4)*ff4_
        odeFunction(t_inner, yy_inner, odeParams, ff5_);

        // Evaluate time derivatives at stage 6:
        t_inner = ti + (a[5] * h);
        _dynorb_rvcopy(sys_size, yyi_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 5, 0)), ff1_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 5, 1)), ff2_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 5, 2)), ff3_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 5, 3)), ff4_, yy_inner);
        _dynorb_raxpy(sys_size, (h * _dynorb_rel(b, 6, 5, 4)), ff4_, yy_inner); // yy_inner = yyi_ + b(6,1)*ff1_ + b(6,2)*ff2_ + b(6,3)*ff3_ + b(6,4)*ff4_ + b(6,5)*ff5_
        odeFunction(t_inner, yy_inner, odeParams, ff6_);

        // Truncation Error Vector:
        _dynorb_rvcopy(sys_size, ff1_, ee);
        _dynorb_rscal(sys_size, (c[0] - c_star[0]), ee);
        _dynorb_raxpy(sys_size, (c[1] - c_star[1]), ff2_, ee);
        _dynorb_raxpy(sys_size, (c[2] - c_star[2]), ff3_, ee);
        _dynorb_raxpy(sys_size, (c[3] - c_star[3]), ff4_, ee);
        _dynorb_raxpy(sys_size, (c[4] - c_star[4]), ff5_, ee);
        _dynorb_raxpy(sys_size, (c[5] - c_star[5]), ff6_, ee);
        _dynorb_rscal(sys_size, h, ee);

        // Maximum Truncation Error:
        trunc_err_max = _dynorb_max_abs_component(sys_size, ee);
        // Allowable Truncation Error:
        trunc_err_allowed = _dynorb_max_abs_component(sys_size, yy);
        trunc_err_allowed = tol * _dynorb_rmax(trunc_err_allowed, 1.0);
        // Fractional Change in step-size:
        delta = (real)pow((trunc_err_allowed / (trunc_err_max + DBL_EPSILON)), (1 / 5)); // COULD GIVE PRROBLEMS float/double

        // If truncation error is within bound --> update solution:
        if (trunc_err_max <= trunc_err_allowed)
        {
            h = _dynorb_rmin(h, (t1 - t));
            t += h;
            // Compute updated solution
            _dynorb_rvcopy(sys_size, ff1_, yy);
            _dynorb_rscal(sys_size, c[0], yy);
            _dynorb_raxpy(sys_size, c[1], ff2_, yy);
            _dynorb_raxpy(sys_size, c[2], ff3_, yy);
            _dynorb_raxpy(sys_size, c[3], ff4_, yy);
            _dynorb_raxpy(sys_size, c[4], ff5_, yy);
            _dynorb_raxpy(sys_size, c[5], ff6_, yy);
            _dynorb_rscal(sys_size, h, yy);
            _dynorb_raxpy(sys_size, 1.0, yyi_, yy);
            // Store the results into solution arrays / realloc if needed:
            if (step >= n_steps)
            {
                // Double the current number of steps (increase size)
                // n_steps *= 2;
                ++n_steps;
                // Update Structure:
                solver_configuration->n_steps = n_steps;

                // Reallocate memory for tt (time array)
                tt = (double *)realloc(tt, n_steps * sizeof(double));
                if (tt == NULL)
                {
                    fprintf(stderr, "Error reallocating memory for tt\n");
                    exit(EXIT_FAILURE);
                }

                // Reallocate memory for YY_t (solution array)
                YY_t = (double *)realloc(YY_t, n_steps * sys_size * sizeof(double));
                if (YY_t == NULL)
                {
                    fprintf(stderr, "Error reallocating memory for YY_t\n");
                    exit(EXIT_FAILURE);
                }
            }
            // Store the results
            tt[step] = t;                                         // Store the current time
            _dynorb_rvcopy(sys_size, yy, &YY_t[step * sys_size]); // Store the current state

            // Increase step:
            ++step;
            // printf("End time step t = %f | h = %f\n", t, h);
        }
        else
        {
            printf("Maximum truncation error = %f | Allowed Truncation error = %f\n", trunc_err_max, trunc_err_allowed);
        }
        // Update time step:
        // h = _dynorb_rmin(delta * h, 4.0 * h);
        h = (h / 1.1) + 0.0 * delta;
        if (h < hmin)
        {
            printf("\n\n Warning: Step size (h = %.10f) fell below\n", h);
            printf("its minimum allowable value (hmin = %.10f) at time %f.\n\n", hmin, t);
            break;
        }
    }
}
