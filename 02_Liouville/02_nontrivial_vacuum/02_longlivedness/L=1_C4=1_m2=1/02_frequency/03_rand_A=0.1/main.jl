using DifferentialEquations # for the actual time evolution
using OrdinaryDiffEq # for ODEs
using Plots # for plotting
using Base.Threads # for parallelization
using StaticArrays # somehow needed to use multiple variables in DifferentialEquations.jl

using Plots, LaTeXStrings, Colors
using Plots.PlotMeasures
using LinearAlgebra

using Random, Distributions

using FFTW # discrete Fourier transform

using JLD2 # for file saving

level = "../"

include(joinpath(level, "src/4th-order-FD-stencils.jl"));
include(joinpath(level, "src/evolution_Liouville_larger_cutoff.jl"));
include(joinpath(level, "src/hamiltonian_Liouville.jl"));
include(joinpath(level, "src/initial_data_waves.jl"));
include(joinpath(level, "src/visualisation.jl"));

#   evolution
#   –––––––––

function artisan_evolution_at_resolution(Nx, stableRandomSeed, pModel, pInit, target_time)
    # unpack model parameters
    (mphi2, mchi2, c4, c, epsDiss) = pModel;
    
    # set the spatial discretization
    NboundaryPadding = 2;#Int(div(Nx,2));  # Number of boundary padding points
    dx = 1/(Nx);  # Grid spacing
        
    pGrid = (dx, Nx, NboundaryPadding);
    # reset a combined set of parameters (residual from old structure ... could be modified)
    # TODO: modify to pGrid, pModel, pInit
    p = (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding);

    # set the time span
    tspan = (0, target_time);

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
        saveat = tspan[end]/10^3, #exp.(range(log(tspan[1]), log(tspan[end]), length=10^4))
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

function evolution_at_param(param, current_target_time, current_res_log2, stableRandomSeed)
    
    #stableRandomSeed = 42
    print("persistent random seed: ", stableRandomSeed, "\n")
    
    print("current characteristic frequency: ", param, "\n")
    
    # set monitoring flags
    convergence_maintained = false;
    lower_bound_only = true;
    
    # parameters of the model
    mphi2 = 1;
    mchi2 = 15;
    c4 = 1;
    c = -1.;
    epsDiss = 0;

    # parameters of the initial data
    a0phi = 0; # effectively sets the relative amplitude to the stochastic ID (since Tkin is kept fixed)
    a0chi = a0phi; # effectively sets the relative amplitude to the stochastic ID (since Tkin is kept fixed)
    k0phi = 1;
    k0chi = 2 * k0phi;
    x0phi = 0;
    x0chi = 1/3;

    offsetphi = sqrt(((mphi2 - mchi2)^2 - 8 * mphi2 * c4) / (32 * c4^2))
    offsetchi = - sqrt(((mphi2 - mchi2)^2 - 8 * mchi2 * c4) / (32 * c4^2))

    aStochastic = 0.1;
    mink = param;
    maxk = param + 4;
    
    desiredTkinPhi = NaN;
    desiredTkinChi = NaN;
    
    # set the combined set of parameters 
    pModel = (mphi2, mchi2, c4, c, epsDiss);
    pInit = (
        a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, 
        offsetphi, offsetchi, 
        aStochastic, mink, maxk, stableRandomSeed,
        desiredTkinPhi, desiredTkinChi
    );
     
    # set some tables to store intermediate output
    resTab = [2^i for i in current_res_log2-2:current_res_log2]
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
        print("current resolution: ", res, "\n")
        # run the evolution
        @time (p, sol, hamiltonian, hamPhi, hamChi) = artisan_evolution_at_resolution(
            res, stableRandomSeed, pModel, pInit, current_target_time
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
    # CONVERGENCE
    #############################
    
    dir_path = string("plots/",stableRandomSeed,"/",param)

    # create the directory if it does not yet exist
    if !isdir(dir_path)
        print("Output plot directory does not exist. Creating it ...\n")
        mkpath(dir_path)
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
    
    # determine the index of convergence loss
    loss_of_convergence_index = findfirst(t -> t > loss_of_convergence_time, solTab[end].t)
    if loss_of_convergence_index === nothing
        loss_of_convergence_index = length(solTab[end].t)
    end
    
    #print(10 * hamPhiTab[end][1],"\n")
    #print(hamPhiTab[end][2:10],"\n")
    
    # determine the onset time of the runaway (10-fold increase in either kinetic energy)
    runaway_index_phi = findfirst(
        energy -> abs(energy) > exp(1) * maximum(abs.(hamPhiTab[end][1:Int(div(length(hamPhiTab[end]),4))+1])), 
        hamPhiTab[end]
    )
    if runaway_index_phi === nothing
        runaway_index_phi = length(hamPhiTab[end])
    end
    runaway_index_chi = findfirst(
        energy -> abs(energy) > exp(1) * maximum(abs.(hamChiTab[end][1:Int(div(length(hamChiTab[end]),4))+1])), 
        hamChiTab[end]
    )
    if runaway_index_chi === nothing
        runaway_index_chi = length(hamChiTab[end])
    end
    runaway_index = min(runaway_index_phi, runaway_index_chi)
    runaway_time = solTab[end].t[runaway_index]
    if runaway_time >= solTab[end].t[end]
        print("No runaway detected.\n")
    else
        print("Runaway detected at time t=",runaway_time,"\n")
    end   
    
    #############################
    # GENERATE REMAINING PLOTS IF DESIRED
    #############################
    
    dir_path = string("plots/",stableRandomSeed,"/",param)

    # create the directory if it does not yet exist
    if !isdir(dir_path)
        #print("Directory does not exist. Creating it...")
        mkpath(dir_path)
    end
        
    # plot energy components
    save_energies_plot(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        #loss_of_convergence_time=loss_of_convergence_time
    )
    save_normalised_energies_plot(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        #loss_of_convergence_time=loss_of_convergence_time
    )
    save_difference_in_energies_plot(
        resTab, pTab, solTab, 
        hamiltonianTab, hamPhiTab, hamChiTab, 
        dir_path,
        #loss_of_convergence_time=loss_of_convergence_time
    )
    
    # plot field heatmaps
    save_density_plots(
        solTab[end], pTab[end], pInit,
        dir_path,
        loss_of_convergence_time=loss_of_convergence_time
    )
    
    # save snapshots
    save_snaps(
        solTab[end], pTab[end];
        snap_intervals=Int(round(length(solTab[end])/1)), 
        yrangeVal=1.2,
        dir_path = dir_path
    );
    
    # animate the fields
    save_animation(
        solTab[end][1:max(1,div(loss_of_convergence_index,10^2)):loss_of_convergence_index], 
        pTab[end],
        join([dir_path, "/animation_Nx=", resTab[end], ".gif"])
    );    
#     # animate frequencies
#     save_animation_momentum_space(
#         solTab[end][1:max(1,div(loss_of_convergence_index,10^2)):loss_of_convergence_index], 
#         pTab[end],
#         join([dir_path, "/animation_momentum_space_Nx=", resTab[end], ".gif"])
#     );
    
    print("Finished plotting.", "\n")
    
    #############################
    # SAVE DATA
    #############################
    
    if runaway_time > loss_of_convergence_time
        print("WARNING: Convergence not maintained until onset of runaway. Resolution insufficient.", "\n")
    else
        convergence_maintained = true
        if runaway_time >= solTab[end].t[end]
            print("WARNING: Lower bound only because target time insufficient.", "\n")
        else
            lower_bound_only = false
        end
    end
    
    dir_path = string("dat/",stableRandomSeed)
    if !isdir(dir_path)
        mkpath(dir_path)
    end

    timesteps = solTab[end].t
    stable_until = min(runaway_time, loss_of_convergence_time)

    @save joinpath(pwd(), dir_path, string(param,".jld2")) param stable_until lower_bound_only timesteps hamiltonianTab hamPhiTab hamChiTab

    print("Saved data.", "\n")

    
    return (runaway_time, convergence_maintained, lower_bound_only)
end

#   main()
#   ––––––

#param_base = 1.2
#param_table = reverse([param_base^i for i in -8:8])

param_table = [freq for freq in 1:1:8]

function main()
    
    # some random seed (can be modified at will)
    stableRandomSeed = rand(1:10^7)
    
    # initialise flags
    convergence_maintained = false;
    lower_bound_only = true;
    
    # set abort criteria ...
    highest_res_log2 = 14;
    max_target_time = 2 * 10^4;
    # ... and their initial values
    current_res_log2 = 10;
    current_target_time = 1;
    
    # initialise the runaway time for handover to next param value
    runaway_time = Inf;
    
    # set table of desired param_table (NOTE: links to scaling assumption below)
    param_base = 1
    param_table = [freq for freq in 1:1:8]
    
    # loop over all values in param_table
    for param in param_table
        
        # re-attempt while flags not positive or until abort criteria met
        while (!convergence_maintained||lower_bound_only) && (current_res_log2 <= highest_res_log2) && (current_target_time <= max_target_time)
            # attempt run and obtain flags
            (runaway_time, convergence_maintained, lower_bound_only) = evolution_at_param(
                param, 
                current_target_time, 
                current_res_log2,
                stableRandomSeed
            )
            # update according to obtained flags
            if convergence_maintained                
                if lower_bound_only
                    current_target_time = current_target_time * 4
                    print("Increasing target time to T = ", current_target_time, "\n")
                else
                    current_target_time = min(runaway_time, current_target_time);
                    print("Target time reset to confidently detected runaway time T = ", current_target_time, "\n")
                end
            else
                current_res_log2 = current_res_log2 + 1;
                print("Increasing resolution from N = ", current_res_log2 - 1, " to ", current_res_log2, "\n")
            end
        end
        
        print("PARAM = ", param, " DONE!\n")
        
        # update target time based on the presumed scaling assumption and adapt the target time accordingly
        current_target_time = runaway_time
        print("Updating target time for next param value from T = ", current_target_time, " ... ")
        current_target_time = current_target_time * exp(param_base)     
        print("to T = ", current_target_time, "\n")
        
        # decrease resolution if convergence was maintained in previous step
#         if convergence_maintained
#             current_res_log2 = current_res_log2 - 1;
#             print("Decreasing resolution from N = ", current_res_log2 + 1, " to ", current_res_log2, "\n")
#         end
        
        # check whether it makes sense to go on; otherwise abort 
        if convergence_maintained && lower_bound_only && current_target_time >= max_target_time
            print("ABORT: maximum target time approached in converged simulation; no use to proceed")
            return
        end
        
        # ensure that current params don't exceed the abort criteria for the next step
        current_res_log2 = min(current_res_log2, highest_res_log2)
        current_target_time = min(current_target_time, max_target_time)
        
        # reset the flags
        convergence_maintained = false;
        lower_bound_only = true;
    end

end

main()

#   export .jl for production run
#   –––––––––––––––––––––––––––––

using NBInclude
nbexport("main.jl", "main.ipynb")