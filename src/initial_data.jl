#   Setup: Prepare initial conditions
#   –––––––––––––––––––––––––––––––––


function add_plane_wave_initial_data!(x, y0, y1, y2, y3, p, pInit)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi       
        ) = pInit
    
    k0phi = k0phi * 2*pi/((Nx)*dx);
    omega0phi = sqrt(k0phi^2 + mphi2);
    
    k0chi = k0chi * 2*pi/((Nx)*dx);
    omega0chi = - sqrt(k0chi^2 + mchi2);
    
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the wave
        y0[n] += abs(a0phi) * sin(k0phi * (x[idx] - x0phi))
        y1[n] += omega0phi * a0phi * cos(k0phi * (x[idx] - x0phi))
        y2[n] += abs(a0chi) * sin(k0chi * (x[idx] - x0chi))
        y3[n] += omega0chi * a0chi * cos(k0chi * (x[idx] - x0chi))
    end
    
    return hcat(y0, y1, y2, y3)

end



function add_wavepacket_initial_data!(x, y0, y1, y2, y3, p, pInit)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi       
        ) = pInit
    
    k0phi = k0phi * 2*pi/(Nx*dx);
    omega0phi = sqrt(k0phi^2);
    
    k0chi = k0chi * 2*pi/(Nx*dx);
    omega0chi = sqrt(k0chi^2);
    
    aphi = a0phi
    achi = a0chi
    
    w0phi = Nx*dx/2^5
    w0chi = Nx*dx/2^5
    
    cphi = 1
    cchi = -1
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the Gaussian wave packet
        y0[n] += aphi * exp(- (x[idx] - x0phi)^2 / (2 * w0phi^2)) * cos(k0phi * (x[idx] - x0phi))
        y1[n] += aphi * cphi * (
            (x[idx] - x0phi) / w0phi^2 * cos(k0phi * (x[idx] - x0phi))
            + k0phi * sin(k0phi * (x[idx] - x0phi))
        ) * exp(- (x[idx] - x0phi)^2 / (2 * w0phi^2))
        y2[n] += achi * exp(- (x[idx] - x0chi)^2 / (2 * w0chi^2)) * cos(k0chi * (x[idx] - x0chi))
        y3[n] += achi * cchi * (
            (x[idx] - x0chi) / w0chi^2 * cos(k0chi * (x[idx] - x0chi))
            + k0chi * sin(k0chi * (x[idx] - x0chi))
        ) * exp(- (x[idx] - x0chi)^2 / (2 * w0chi^2))
    end
    
    return hcat(y0, y1, y2, y3)

end


function add_moving_Gaussian_initial_data!(x, y0, y1, y2, y3, p, pInit)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi       
        ) = pInit
    
    Lbox = (Nx*dx)
    
    k0phi = k0phi * 2*pi/Lbox;
    k0chi = k0chi * 2*pi/Lbox;
    
    aphi = a0phi
    achi = a0chi
    
    w0phi = 1/k0phi
    w0chi = 1/k0chi
    
    cphi = 1
    cchi = -1
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the Gaussian wave packet
        y0[n] += aphi * exp(- (x[idx] - x0phi)^2 / (2 * w0phi^2))
        y1[n] += aphi * cphi * (x[idx] - x0phi) / w0phi^2 * exp(- (x[idx] - x0phi)^2 / (2 * w0phi^2))
        y2[n] += achi * exp(- (x[idx] - x0chi)^2 / (2 * w0chi^2))
        y3[n] += achi * cchi * (x[idx] - x0chi) / w0chi^2 * exp(- (x[idx] - x0chi)^2 / (2 * w0chi^2))
    end
    
    return hcat(y0, y1, y2, y3)

end



function add_gaussian_bump_initial_data!(x, y0, y1, y2, y3, p, pInit)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi       
        ) = pInit
    
    k0phi = k0phi * 2*pi/(Nx*dx);
    omega0phi = sqrt(k0phi^2);
    
    k0chi = k0chi * 2*pi/(Nx*dx);
    omega0chi = sqrt(k0chi^2);
    
    aphi = a0phi
    achi = a0chi
    
    w0phi = Nx*dx/2^5
    w0chi = Nx*dx/2^5
    
    cphi = 1
    cchi = -1
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the Gaussian wave packet
        y0[n] += aphi * exp(- (x[idx] - x0phi)^2 / (2 * w0phi^2)) * cos(k0phi * (x[idx] - x0phi))
        y2[n] += achi * exp(- (x[idx] - x0chi)^2 / (2 * w0chi^2)) * cos(k0chi * (x[idx] - x0chi))
    end
    
    return hcat(y0, y1, y2, y3)

end




function add_offset!(x, y0, y1, y2, y3, p, pInit)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi      
        ) = pInit
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        
        # add the offset
        
        y0[n] += offsetphi
        
        y2[n] += offsetchi
    
    end
    
    return hcat(y0, y1, y2, y3)

end


function add_noise!(x, y0, y1, y2, y3, p, noise_amplitude)
    
    # unpack
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        
        # add the offset
        
        y0[n] += noise_amplitude * rand()
        y1[n] += noise_amplitude * rand()        
        y2[n] += noise_amplitude * rand()
        y3[n] += noise_amplitude * rand()
    
    end
    
    return hcat(y0, y1, y2, y3)

end




function add_stochastic_initial_data_at_wavenumber!(x, y0, y1, y2, y3, p, pInit, k)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi     
        ) = pInit
    
    Lbox= ((Nx)*dx)
    
    kphi = k*2*pi/Lbox
    kchi = k*2*pi/Lbox
    omegaphi = sqrt(kphi^2 + mphi2)
    omegachi = sqrt(kchi^2 + mchi2)

    sigmaphik = sqrt(1/(2*kphi*Lbox))
    sigmachik = sqrt(1/(2*kchi*Lbox))

    deltaphikR = rand(Normal(0,sigmaphik))
    deltachikR = rand(Normal(0,sigmachik))
    deltaphikL = rand(Normal(0,sigmaphik))
    deltachikL = rand(Normal(0,sigmachik))
    
    x0phi = Lbox*rand()
    x0chi = Lbox*rand()
    x0phi = Lbox*rand()
    x0chi = Lbox*rand()
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the right-moving one
        y0[n] += aStochastic * deltaphikR * sin(kphi * (x[idx] - x0phi))        
        y1[n] += aStochastic * deltaphikR * omegaphi * cos(kphi * (x[idx] - x0phi))
        y2[n] += aStochastic * deltachikR * sin(kchi * (x[idx] - x0chi))
        y3[n] += aStochastic * deltachikR * omegachi * cos(kchi * (x[idx] - x0chi))
        # add the left-moving one
        y0[n] += aStochastic * deltaphikL * sin(kphi * (x[idx] - x0phi))        
        y1[n] += - aStochastic * deltaphikL * omegaphi * cos(kphi * (x[idx] - x0phi))
        y2[n] += aStochastic * deltachikL * sin(kchi * (x[idx] - x0chi))
        y3[n] += - aStochastic * deltachikL * omegachi * cos(kchi * (x[idx] - x0chi))
    end

end



function add_stochastic_initial_data_at_wavenumber_phi_only!(x, y0, y1, y2, y3, p, pInit, k)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi     
        ) = pInit
    
    Lbox= ((Nx)*dx)
    
    kphi = k*2*pi/Lbox
    omegaphi = sqrt(kphi^2 + mphi2)

    sigmaphik = sqrt(1/(2*kphi*Lbox))

    deltaphikR = rand(Normal(0,sigmaphik))
    deltaphikL = rand(Normal(0,sigmaphik))
    
    x0phi = Lbox*rand()
    x0chi = Lbox*rand()
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the right-moving one
        y0[n] += aStochastic * deltaphikR * sin(kphi * (x[idx] - x0phi))        
        y1[n] += aStochastic * deltaphikR * omegaphi * cos(kphi * (x[idx] - x0phi))
        # add the left-moving one
        y0[n] += aStochastic * deltaphikL * sin(kphi * (x[idx] - x0phi))        
        y1[n] += - aStochastic * deltaphikL * omegaphi * cos(kphi * (x[idx] - x0phi))
    end

end


function add_stochastic_initial_data_at_wavenumber_chi_only!(x, y0, y1, y2, y3, p, pInit, k)
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi     
        ) = pInit
    
    Lbox= ((Nx)*dx)
    
    kchi = k*2*pi/Lbox
    omegachi = sqrt(kchi^2 + mchi2)

    sigmachik = sqrt(1/(2*kchi*Lbox))

    deltachikR = rand(Normal(0,sigmachik))
    deltachikL = rand(Normal(0,sigmachik))
    
    x0phi = Lbox*rand()
    x0chi = Lbox*rand()
        
    for n = NboundaryPadding + 1:Nx + NboundaryPadding
        # set the right index in the x-range (account for boundary)
        idx = n - NboundaryPadding
        # add the right-moving one
        y2[n] += aStochastic * deltachikR * sin(kchi * (x[idx] - x0chi))
        y3[n] += aStochastic * deltachikR * omegachi * cos(kchi * (x[idx] - x0chi))
        # add the left-moving one
        y2[n] += aStochastic * deltachikL * sin(kchi * (x[idx] - x0chi))
        y3[n] += - aStochastic * deltachikL * omegachi * cos(kchi * (x[idx] - x0chi))
    end

end
