using DifferentialEquations # for the actual time evolution
using OrdinaryDiffEq # for ODEs
using Plots # for plotting
using FFMPEG
using Base.Threads # for parallelization
using StaticArrays # somehow needed to use multiple variables in DifferentialEquations.jl

using Plots, LaTeXStrings, Colors
using Plots.PlotMeasures
using LinearAlgebra

using Random, Distributions

using FFTW # discrete Fourier transform

level = "../../../../../"

include(joinpath(level, "src/4th-order-FD-stencils.jl"));
include(joinpath(level, "src/evolution_simple_lambda22.jl"));
include(joinpath(level, "src/hamiltonian_simple_lambda22.jl"));
include(joinpath(level, "src/initial_data_waves.jl"));
include(joinpath(level, "src/visualisation.jl"));

# modified plotting routines for the paper
include("../../visualisation_paper-version.jl");

#   main: evolution
#   –––––––––––––––

function artisan_evolution_at_resolution(Nx, stableRandomSeed, pModel, pInit)
    # unpack model parameters
    (mphi2, mchi2, lambdaCross, lambdaSelf, sigma) = pModel;
    
    # set the spatial discretization
    NboundaryPadding = 2;  # Number of boundary padding points
    dx = 1/Nx;  # Grid spacing
    
    pGrid = (dx, Nx, NboundaryPadding);
    # reset a combined set of parameters (residual from old structure ... could be modified)
    # TODO: modify to pGrid, pModel, pInit
    p = (dx, mphi2, mchi2, lambdaCross, lambdaSelf, sigma, Nx, NboundaryPadding);

    # set the time span
    tspan = (0, 2 * 1.1);

    # set the evolution method
    time_integration_method = RK4();
    
    # generate initial conditions
    u0 = initial_data(
        range(0, step=dx, length=(Nx + 2 * NboundaryPadding)),
        p, 
        pInit
    );
        
    # set the problem
    prob = ODEProblem(finite_differenced_pde_with_bc!, u0, tspan, p);

    sol = solve(
        prob, time_integration_method, 
        saveat = tspan[end]/10^3,
#         saveat = vcat(
#             exp.(range(log(0.001), log(1), length=10^3)),
#             range(tspan[1], tspan[end], length=10^3)
#         ),
        dt=dx/4, 
        adaptive = false, 
        dense=false, 
        maxiters=typemax(Int),
        callback=field_size_callback
    );
        
    # obtain the hamiltonian
    hamiltonian = zeros(length(sol.u))
    hamphi = zeros(length(sol.u))
    hamchi = zeros(length(sol.u))
    for i = 1:length(sol.u)
        hamiltonian[i] = nintegrate_simps(hamiltonian_density(sol.u[i], p), dx)
        hamphi[i] = nintegrate_simps(hamiltonian_phi(sol.u[i], p), dx)
        hamchi[i] = nintegrate_simps(hamiltonian_chi(sol.u[i], p), dx)
    end
    
    return (p, sol, hamiltonian, hamphi, hamchi)
end

function main()
    # parameters of the model
    mphi2 = 0;
    mchi2 = 10000;
    lambdaCross = + 100.0;
    lambdaSelf = 0.0;
    sigma = - 1;

    # parameters of the initial data
    a0phi = + 1.0; # effectively sets the relative amplitude to the stochastic ID (since Tkin is kept fixed)
    a0chi = a0phi; # effectively sets the relative amplitude to the stochastic ID (since Tkin is kept fixed)
    k0phi = 1;
    k0chi = 5 * k0phi;
    x0phi = 0;
    x0chi = 1/3;

    offsetphi = 0;
    offsetchi = 0;

    aStochastic = 0; # effectively sets the relative amplitude to the plane wave ID
    mink = 1;
    maxk = 4;
    
    desiredValuePhi = NaN; # no rescaling desired
    desiredValueChi = NaN; # no rescaling desired
    
    # some random seed (can be modified at will)
    stableRandomSeed = 42#rand(1:10^7)
    #stableRandomSeed = 313371
    print("persistent random seed: ", stableRandomSeed, "\n")

    # set the combined set of parameters 
    pModel = (mphi2, mchi2, lambdaCross, lambdaSelf, sigma);
    pInit = (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, 
        offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredValuePhi, desiredValueChi
    );
     
    # set some tables to store output
    resTab = [2^i for i in 7:11]
    pTab = []
    solTab = []
    hamiltonianTab = []
    hamPhiTab = []
    hamChiTab = []
    
    #############################
    # evolution
    #############################

    # run evolution
    for res in resTab
        print("current resolution: ", res, " ... \n")
        # run the evolution
        @time (p, sol, hamiltonian, hamPhi, hamChi) = artisan_evolution_at_resolution(
            res, stableRandomSeed, pModel, pInit
        )
        print("... terminated", "\n")
        
        # unpack parameters
        (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
        # append the results
        push!(pTab, p)
        push!(solTab, sol)
        push!(hamiltonianTab, hamiltonian)
        push!(hamPhiTab, hamPhi)
        push!(hamChiTab, hamChi)
    end
    
    
    #############################
    # PLOTTING
    #############################
    
    
    dir_path = "plots"

    # create the directory if it does not yet exist
    if !isdir(dir_path)
        print("Output plot directory does not exist. Creating it ...\n")
        mkdir(dir_path)
    else
        print("Output plot directory already exists.\n")
    end
    
    # plot and determine convergence 
    loss_of_convergence_time = save_convergence_plots(
        resTab, pTab, solTab, 
        hamiltonianTab, 
        dir_path
    )
    if loss_of_convergence_time >= solTab[end].t[end]
        print("Convergence kept at all times.\n")
    else
        print("Convergence lost at time t=",loss_of_convergence_time,"\n")
    end
    # loss_of_convergence_time = solTab[end].t[end]
    
    # plot energy components
    save_energies_plot(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )
#     save_normalised_energies_plot(
#         resTab, pTab, solTab, 
#         hamiltonianTab, hamPhiTab, hamChiTab, 
#         dir_path,
#         loss_of_convergence_time=loss_of_convergence_time
#     )
    save_difference_in_energies_plot(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )    
    save_energies_plot_paper(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )
    save_difference_in_energies_plot_paper(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )  
    
    # plot field heatmaps
    save_density_plots(
        solTab[end], pTab[end], pInit,
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )    
    save_density_plots_paper(
        solTab[end], pTab[end], pInit,
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )
    
    save_snaps(
        solTab[end], pTab[end];
        snap_intervals=Int(round(length(solTab[end])/1)), 
        yrangeVal=1.2,
        dir_path = dir_path
    );
    
    # determine the index of convergence loss
    loss_of_convergence_index = findfirst(t -> t > loss_of_convergence_time, solTab[end].t)
    if loss_of_convergence_index === nothing
        loss_of_convergence_index = length(solTab[end].t)
    end
    # and then animate the fields
    save_animation(
        solTab[end],#[1:max(1,div(loss_of_convergence_index,10^2)):loss_of_convergence_index], 
        pTab[end],
        dir_path,
        fps=16
    );
    # and the frequencies
#     save_animation_momentum_space(
#         solTab[end][1:max(1,div(loss_of_convergence_index,10^3)):loss_of_convergence_index], 
#         pTab[end],
#         dir_path
#     );
    
    print("Run finished.")
end

main()

#   output
#   ––––––

using NBInclude
nbexport("main.jl", "main.ipynb")