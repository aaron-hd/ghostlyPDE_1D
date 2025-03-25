#   setup: calculate the integrated Hamiltonian
#   –––––––––––––––––––––––––––––––––––––––––––

#   Here is a direct implementation of Simpson's rule. (adapted from
#   https://mmas.github.io/simpson-integration-julia ).

using NumericalIntegration

#   Now define a function that determines the Hamiltonian density.

function hamiltonian_density(fields_at_t, p)
    # unpack parameters
    (dx, mphi2, mchi2, lambdaGhost, lambdaSelf, sigma, Nx, NboundaryPadding) = p
    
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
            + sigma * 1//2 * y3[idx]^2
            + sigma * 1//2 * first_derivative(y2, idx, dx)^2
            # masses
            + mphi2 * y0[idx]^2 / 2 
            + sigma * mchi2 * y2[idx]^2 / 2 
            # interactions
            + lambdaGhost * y0[idx]^2 * y2[idx]^2
            + lambdaSelf * y0[idx]^4 
            + sigma * lambdaSelf * y2[idx]^4
        )
    end
    ham
end

function hamiltonian_phi(fields_at_t, p)
    # unpack parameters
    (dx, mphi2, mchi2, lambdaGhost, lambdaSelf, sigma, Nx, NboundaryPadding) = p
    
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
            # mass phi
            + mphi2 * y0[idx]^2 / 2 
            # self-interaction
            + lambdaSelf * y0[idx]^4
        )
    end
    hamphi
end

function hamiltonian_chi(fields_at_t, p)
    # unpack parameters
    (dx, mphi2, mchi2, lambdaGhost, lambdaSelf, sigma, Nx, NboundaryPadding) = p
    
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
            + sigma * 1//2 * y3[idx]^2
            + sigma * 1//2 * first_derivative(y2, idx, dx)^2
            # mass phi
            + sigma * mchi2 * y2[idx]^2 / 2 
            # self-interactions
            + sigma * lambdaSelf * y2[idx]^4
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
