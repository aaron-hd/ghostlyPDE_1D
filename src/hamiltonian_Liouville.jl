#   setup: calculate the integrated Hamiltonian
#   –––––––––––––––––––––––––––––––––––––––––––

#   Here is a direct implementation of Simpson's rule. (adapted from
#   https://mmas.github.io/simpson-integration-julia ).

using NumericalIntegration

#   Now define a function that determines the Hamiltonian density.

function hamiltonian_density(fields_at_t, p)
    # unpack parameters
    dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding = p
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y0, y1, y2, y3 = [
            fields_at_t[1:N], 
            fields_at_t[N+1:2N], 
            fields_at_t[2N+1:3N], 
            fields_at_t[3N+1:4N]
        ]
    
    # prepare and calculate the hamiltonian density    
    # NOTE: to achieve periodicity for the integrator, 
    # we need to include one of the boundary points
    ham = zeros(Nx + 1)  # Include the extra point for periodicity
    for n = NboundaryPadding + 1:Nx + NboundaryPadding + 1
        idx = mod(n - 1, Nx) + NboundaryPadding + 1
        ham[n - NboundaryPadding] = (
            # phi
            + 1//2 * y1[idx]^2
            + 1//2 * first_derivative(y0, idx, dx)^2
            # chi
            - 1//2 * y3[idx]^2
            - 1//2 * first_derivative(y2, idx, dx)^2
            # cross-coupling (including masses)
            + (mphi2*y0[idx]^2 
            - mchi2*y2[idx]^2)/2 
            - ((-mchi2 + mphi2)/(2*c) * (y0[idx]^2 - y2[idx]^2)^2)
            + c4*(y0[idx]^2 - y2[idx]^2)^3 
            + c*c4*(-y0[idx]^4 + y2[idx]^4)
        )
    end
    ham
end

function hamiltonian_phi(fields_at_t, p)
    # unpack parameters
    dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding = p
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y0, y1 = [
            fields_at_t[1:N], 
            fields_at_t[N+1:2N]
        ]
    
    # prepare and calculate the hamiltonian density
    hamphi = zeros(Nx + 1)
    for n = NboundaryPadding + 1:Nx + NboundaryPadding + 1
        idx = mod(n - 1, Nx) + NboundaryPadding + 1
        hamphi[n - NboundaryPadding] = (
            # phi
            + 1//2 * y1[idx]^2
            + 1//2 * first_derivative(y0, idx, dx)^2
            # mass
            + (mphi2*y0[idx]^2)/2
            # self-interactions
            - (-mchi2 + mphi2)/(2*c) * y0[idx]^4
            + c4*y0[idx]^6 
            - c*c4*y0[idx]^4
        )
    end
    hamphi
end

function hamiltonian_chi(fields_at_t, p)
    # unpack parameters
    dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding = p
    
    # unpack fields
    N = size(fields_at_t, 1)    
    y2, y3 = [
            fields_at_t[2N+1:3N], 
            fields_at_t[3N+1:4N]
        ]
    
    # prepare and calculate the hamiltonian density
    hamchi = zeros(Nx + 1)
    for n = NboundaryPadding + 1:Nx + NboundaryPadding + 1
        idx = mod(n - 1, Nx) + NboundaryPadding + 1
        hamchi[n-NboundaryPadding] = (
            # kinetic chi
            - 1//2 * y3[idx]^2
            - 1//8 * first_derivative(y2, idx, dx)^2
            # mass chi
            - (mchi2*y2[idx]^2)/2 
            # self-interactions
            - (-mchi2 + mphi2)/(2*c) * y2[idx]^4
            - c4*y2[idx]^6
            + c*c4*y2[idx]^4
        )
    end
    hamchi
end



function nintegrate_simps(y::Vector, h::Number)
    n = length(y)-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    s = - sum(y[1:2:n] + 4*y[2:2:n] + y[3:2:n+1])
    return h/3 * s
end;

function nintegrate_simps38(y::Vector, h::Number)
    n = length(y) - 1  # Number of intervals
    
    if n % 3 == 0
        # Apply Simpson's 3/8 Rule across the entire domain
        s = y[1] + y[end]  # First and last terms
        
        # Add terms grouped by 3/8 Rule weights
        s += 3 * sum(y[2:3:n] .+ y[3:3:n])  # Terms with weight 3
        s += 2 * sum(y[4:3:n-1])            # Terms with weight 2

        return (3 * h / 8) * s
    else
        # Hybrid approach: Apply Simpson's 3/8 Rule for largest multiple of 3 intervals
        m = n - (n % 3)  # Largest multiple of 3 less than or equal to n

        # 3/8 Rule over the first m intervals
        s = y[1] + y[m+1]  # Include first and (m+1)-th term
        s += 3 * sum(y[2:3:m] .+ y[3:3:m])  # Terms with weight 3
        s += 2 * sum(y[4:3:m-1])            # Terms with weight 2

        simps38_part = (3 * h / 8) * s

        # Remaining intervals handled by Trapezoidal Rule
        remainder_part = 0.0
        if m < n
            for i in m+1:n
                remainder_part += (y[i] + y[i+1]) * (h / 2)
            end
        end

        return simps38_part + remainder_part
    end
end
