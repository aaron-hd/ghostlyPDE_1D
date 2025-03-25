#   Setup: Prepare the ODE system (i.e., the PDE after finite differencing
#   in spatial coordinate)
#   ––––––––––––––––––––––––––

function finite_differenced_pde_with_bc!(du, u, p, t)
    # Unpack u into separate arrays
    y0, y1, y2, y3 = [u[:, i] for i in 1:4]
    
    # Unpack du into separate arrays
    dy0dt, dy1dt, dy2dt, dy3dt = [du[:, i] for i in 1:4]
    
    dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding = p
    
    # update non-boundary points
    Threads.@threads for n = NboundaryPadding + 1:Nx + NboundaryPadding
        dy0dt[n] = (
            # kinetic term
            y1[n]
        )

        dy1dt[n] = (
            # kinetic term
            second_derivative(y0, n, dx)
            # potential terms
            - mphi2 * y0[n] 
            + 4 * c * c4 * y0[n]^3 
            + (
                2 / c * (mphi2 - mchi2) * y0[n] * (
                    + y0[n]^2 
                    - y2[n]^2
                )
            ) 
            - 6 * c4 * y0[n] * (y0[n]^2 - y2[n]^2)^2
        )

        dy2dt[n] = (
            # kinetic term
            y3[n]
        )

        dy3dt[n] = (
            # kinetic term
            second_derivative(y2, n, dx)
            # potential terms
            - mchi2 * y2[n] 
            + 4 * c * c4 * y2[n]^3 
            + (
                2 / c * (mphi2 - mchi2) * y2[n] * (
                    + y0[n]^2 
                    - y2[n]^2
                )
            ) 
            - 6 * c4 * y2[n] * (y0[n]^2 - y2[n]^2)^2
        )
    end
    
    # update boundary points (periodic boundary conditions)
    for n in 1:NboundaryPadding
        dy0dt[n] = dy0dt[Nx + n]
        dy1dt[n] = dy1dt[Nx + n]
        dy2dt[n] = dy2dt[Nx + n]
        dy3dt[n] = dy3dt[Nx + n]

        dy0dt[Nx + NboundaryPadding + n] = dy0dt[NboundaryPadding + n]
        dy1dt[Nx + NboundaryPadding + n] = dy1dt[NboundaryPadding + n]
        dy2dt[Nx + NboundaryPadding + n] = dy2dt[NboundaryPadding + n]
        dy3dt[Nx + NboundaryPadding + n] = dy3dt[NboundaryPadding + n]
    end
    
    # Concatenate arrays again and update
    du .= hcat(dy0dt, dy1dt, dy2dt, dy3dt)
end

# Define a callback that checks if any field in 'u' grows too large

function check_field_too_large(u, t, integrator)
    max_value_allowed = 10^3  # Define your threshold here
    return any(x -> x > max_value_allowed, u)  # True if any field exceeds threshold
end

function terminate_if_large!(integrator)
    println("Terminating because one of the fields grew too large at time t = $(integrator.t).")
    terminate!(integrator)
end

# Create a DiscreteCallback that triggers if the check_field_too_large returns true
field_size_callback = DiscreteCallback(
    check_field_too_large, # Condition to check field size
    terminate_if_large!    # Action to terminate the solver
);
