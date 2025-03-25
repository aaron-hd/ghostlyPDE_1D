#   Setup: Prepare initial conditions
#   –––––––––––––––––––––––––––––––––


include("./initial_data.jl")


function initial_data(x, p, pInit)
    
    # unpack
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi        
        ) = pInit
    
    # initialise the field variables
    y0 = zeros((Nx + 2 * NboundaryPadding))
    y1 = zeros((Nx + 2 * NboundaryPadding))
    y2 = zeros((Nx + 2 * NboundaryPadding))
    y3 = zeros((Nx + 2 * NboundaryPadding))
    
    # add the non-stochastic initial data
    add_plane_wave_initial_data!(x, y0, y1, y2, y3, p, pInit)
    
    # add the stochastic initial data provided the suitable persistent seed
    Random.seed!(stableRandomSeed)
    
    for k = mink:maxk
        
        add_stochastic_initial_data_at_wavenumber!(x, y0, y1, y2, y3, p, pInit, k)       
    
    end
    
    # calculate current maximum amplitude
    ampphi = maximum(abs.(y0))
    ampchi = maximum(abs.(y2))
    
    if (desiredValuePhi isa Number && desiredValueChi isa Number) && !(isnan(desiredValuePhi) || isnan(desiredValueChi))
        
        print("\t", "user-assigned rescaling:", "\n")
        print("\t", "\t", "max |amplitude| phi before rescaling: ", ampphi, "\n")
        print("\t", "\t", "max |amplitude| chi before rescaling: ", ampchi, "\n")

        # rescale each component field
        y0 = y0*abs(desiredValuePhi/ampphi)
        y1 = y1*abs(desiredValuePhi/ampphi)
        y2 = y2*abs(desiredValueChi/ampchi)
        y3 = y3*abs(desiredValueChi/ampchi)

        # recalculate current maximum amplitude
        ampphi = maximum(abs.(y0))
        ampchi = maximum(abs.(y2))

        print("\t", "\t", "max |amplitude| phi after rescaling: ", ampphi, "\n")
        print("\t", "\t", "max |amplitude| chi after rescaling: ", ampchi, "\n")
    
    else
    
        print("\t", "user-assigned no rescaling:", "\n")
        print("\t", "\t", "max |amplitude| phi before rescaling: ", ampphi, "\n")
        print("\t", "\t", "max |amplitude| chi before rescaling: ", ampchi, "\n")
    
    end
    
    # add the offset
    add_offset!(x, y0, y1, y2, y3, p, pInit)
    
    # apply periodic boundary conditions
    for n in 1:NboundaryPadding
        y0[n] = y0[Nx + n]
        y1[n] = y1[Nx + n]
        y2[n] = y2[Nx + n]
        y3[n] = y3[Nx + n]

        y0[Nx + NboundaryPadding + n] = y0[NboundaryPadding + n]
        y1[Nx + NboundaryPadding + n] = y1[NboundaryPadding + n]
        y2[Nx + NboundaryPadding + n] = y2[NboundaryPadding + n]
        y3[Nx + NboundaryPadding + n] = y3[NboundaryPadding + n]
    end
    
    return hcat(y0, y1, y2, y3)

end

