phicolor = colorant"#c51b8a";
chicolor = colorant"#00FFFF";






function save_energies_plot_paper(resTab, pTab, solTab, hamiltonianTab, hamPhiTab, hamChiTab, dir_path; loss_of_convergence_time=Inf)
    
    # determine limit values for plots
    if loss_of_convergence_time isa Number
        tlim = (solTab[end].t[2], loss_of_convergence_time)
    else
        tlim = (solTab[end].t[2], solTab[end].t[end])
    end
    loss_of_convergence_time = min(solTab[end].t[end],loss_of_convergence_time)
    dt = solTab[end].t[2] - solTab[end].t[1];
    Nt = Int(min(length(solTab[end].t), div(loss_of_convergence_time,dt)));
    ylimSetAuto = (
        maximum([minimum([
#             minimum(abs.(hamiltonianTab[end][1:Nt])),
            minimum(abs.(hamPhiTab[end][1:Nt])),
            minimum(abs.(hamChiTab[end][1:Nt])),
#             minimum(abs.(hamiltonianTab[end][1:Nt] - hamPhiTab[end][1:Nt] - hamChiTab[end][1:Nt]))     
        ]),10^-6]),
        minimum([maximum([
#             maximum(abs.(hamiltonianTab[end][1:Nt])),
            maximum(abs.(hamPhiTab[end][1:Nt])),
            maximum(abs.(hamChiTab[end][1:Nt])),
#             maximum(abs.(hamiltonianTab[end][1:Nt] - hamPhiTab[end][1:Nt] - hamChiTab[end][1:Nt]))
        ]),10^11])
    );
    
    ## PLOT ONCE WITHOUT ABSOLUTE VALUE
    
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 300, size = (440,600), 
        legend=:topleft,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        #ylims=ylimSetAuto,
        xlabel=L"$t/L$", 
        ylabel=L"$\textrm{kinetic}\;\;\textrm{energies}$",
        left_margin = 20px,
        right_margin = 20px
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #print(dH[end],"\n")
        (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
        # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian,
#             label=latexstring("H^{N=",res,"}"),
#             linewidth=2, linecolor="black", linestyle=:dash
#         )
        plot!(
            energiesPlot, 
            sol.t, hamPhi,
#             label=latexstring("H_\\phi^{N=",res,"}"),
            label=latexstring("H_\\phi"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, hamChi,
#             label=latexstring("H_\\chi^{N=",res,"}"),
            label=latexstring("H_\\chi"),
            linewidth=2, linecolor=chicolor
        )
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian-hamPhi-hamChi,
#             label=latexstring("V_\\mathrm{int}^{N=",res,"}"),
#             linewidth=2, linecolor="orange", linestyle=:dot
#         )
    end

#     for i = 1:(length(resTab)-1)
#         # unpack parameters
#         p = pTab[i]
#         res = resTab[i]
#         sol = solTab[i]
#         hamiltonian = hamiltonianTab[i]
#         hamPhi = hamPhiTab[i]
#         hamChi = hamChiTab[i]
#         #print(dH[end],"\n")
#         (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
#         # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian,
#             label=false,#latexstring("H^{N=",res,"}"),
#             linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dash
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, hamPhi,
#             label=false,#latexstring("H_\\phi^{N=",res,"}"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, hamChi,
#             label=false,#latexstring("H_\\chi^{N=",res,"}"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian-hamPhi-hamChi,
#             label=false,#latexstring("V_\\mathrm{int}^{N=",res,"}"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
#         )
#     end

    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "paper_energies.png"));
    
    ## RE-PLOT WITH ABSOLUTE VALUE FOR LOG-PLOTS
    
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 300, size = (440,600), 
        legend=:topleft,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        ylims=ylimSetAuto,
        xlabel=L"$t/L$", 
        ylabel=L"$\textrm{kinetic}\;\;\textrm{energies}$",
        left_margin = 20px,
        right_margin = 20px
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #print(dH[end],"\n")
        (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
        # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian),
#             label=latexstring("|H^{N=",res,"}|"),
#             linewidth=2, linecolor="black", linestyle=:dash
#         )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamPhi),
#             label=latexstring("|H_\\phi^{N=",res,"}|"),
            label=latexstring("|H_\\phi|"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamChi),
#             label=latexstring("|H_\\chi^{N=",res,"}|"),
            label=latexstring("|H_\\chi|"),
            linewidth=2, linecolor=chicolor
        )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian-hamPhi-hamChi),
#             label=latexstring("|V_\\mathrm{int}^{N=",res,"}|"),
#             linewidth=2, linecolor="orange", linestyle=:dot
#         )
    end

#     for i = 1:(length(resTab)-1)
#         # unpack parameters
#         p = pTab[i]
#         res = resTab[i]
#         sol = solTab[i]
#         hamiltonian = hamiltonianTab[i]
#         hamPhi = hamPhiTab[i]
#         hamChi = hamChiTab[i]
#         #print(dH[end],"\n")
#         (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
#         # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian),
#             label=false,#latexstring("|H^{N=",res,"}|"),
#             linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dash
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamPhi),
#             label=false,#latexstring("|H_\\phi^{N=",res,"}|"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamChi),
#             label=false,#latexstring("|H_\\chi^{N=",res,"}|"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian-hamPhi-hamChi),
#             label=false,#latexstring("|V_\\mathrm{int}^{N=",res,"}|"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
#         )
#     end
    
    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto,
        xlim=tlim
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "paper_energies_log.png"));
    
    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto,
        xscale=:log,
        xlim=tlim
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "paper_energies_loglog.png"));
    
end



function save_difference_in_energies_plot_paper(resTab, pTab, solTab, hamiltonianTab, hamPhiTab, hamChiTab, dir_path; loss_of_convergence_time=Inf)
    
    # determine limit values for plots
    if loss_of_convergence_time isa Number
        tlim = (solTab[end].t[2], loss_of_convergence_time)
    else
        tlim = (solTab[end].t[2], solTab[end].t[end])
    end
    loss_of_convergence_time = min(solTab[end].t[end],loss_of_convergence_time)  
    dt = solTab[end].t[2] - solTab[end].t[1];
    Nt = Int(min(length(solTab[end].t), div(loss_of_convergence_time,dt)));
    ylimSetAuto = (
        maximum([minimum([
#             minimum(abs.(hamiltonianTab[end][2:Nt] .- hamiltonianTab[end][1])),
            minimum(abs.(hamPhiTab[end][2:Nt] .- hamPhiTab[end][1])),
            minimum(abs.(hamChiTab[end][2:Nt] .- hamChiTab[end][1]))    
        ]),10^-6]),
        minimum([maximum([
#             maximum(abs.(hamiltonianTab[end][2:Nt] .- hamiltonianTab[end][1])),
            maximum(abs.(hamPhiTab[end][2:Nt] .- hamPhiTab[end][1])),
            maximum(abs.(hamChiTab[end][2:Nt] .- hamChiTab[end][1]))  
        ]),10^11])
    );
   
    ## PLOT ONCE WITHOUT ABSOLUTE VALUE
    
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 300, size = (440,600), 
        legend=:topleft,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        ylims=ylimSetAuto,
        xlabel=L"$t/L$", 
        ylabel=L"$\textrm{kinetic}\;\;\textrm{energies}$",
        left_margin = 20px,
        right_margin = 20px
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #build difference
        hamiltonian = hamiltonian.-hamiltonian[1]
        hamPhi = hamPhi.-hamPhi[1]
        hamChi = hamChi.-hamChi[1]
        # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian,
#             label=latexstring("H^{N=",res,"} - H^{0}"),
#             linewidth=2, linecolor="black", linestyle=:dash
#         )
        plot!(
            energiesPlot, 
            sol.t, hamPhi,
#             label=latexstring("H_\\phi^{N=",res,"} - H_\\phi^{0}"),
            label=latexstring("H_\\phi - H_\\phi^{0}"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, hamChi,
#             label=latexstring("H_\\chi^{N=",res,"} - H_\\chi^{0}"),
            label=latexstring("H_\\chi - H_\\chi^{0}"),
            linewidth=2, linecolor=chicolor
        )
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian-hamPhi-hamChi,
#             label=latexstring("V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}"),
#             linewidth=2, linecolor="orange", linestyle=:dot
#         )
    end

#     for i = 1:(length(resTab)-1)
#         # unpack parameters
#         p = pTab[i]
#         res = resTab[i]
#         sol = solTab[i]
#         (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
#         hamiltonian = hamiltonianTab[i]
#         hamPhi = hamPhiTab[i]
#         hamChi = hamChiTab[i]
#         #build difference
#         hamiltonian = hamiltonian.-hamiltonian[1]
#         hamPhi = hamPhi.-hamPhi[1]
#         hamChi = hamChi.-hamChi[1]
#         # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian,
#             label=false,#latexstring("H^{N=",res,"} - H^{0}"),
#             linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dash
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, hamPhi,
#             label=false,#latexstring("H_\\phi^{N=",res,"} - H_\\phi^{0}"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, hamChi,
#             label=false,#latexstring("H_\\chi^{N=",res,"} - H_\\chi^{0}"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, hamiltonian-hamPhi-hamChi,
#             label=false,#latexstring("V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
#         )
#     end
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "paper_energies_difference.png"));
    
    ## RE-PLOT WITH ABSOLUTE VALUE FOR LOG-PLOTS
   
    #set up the plot environment
    energiesPlot = plot(
        layout=(1), framestyle=:box, dpi = 300, size = (440,600), 
        legend=:topleft,
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{covergence plot Hamiltonian\\;\\;}H"),
        xlim=tlim,
        ylims=ylimSetAuto,
        xlabel=L"$t/L$", 
        ylabel=L"$\textrm{kinetic}\;\;\textrm{energies}$",
        left_margin = 20px,
        right_margin = 20px
    )

    for i = length(resTab):length(resTab)
        # unpack parameters
        p = pTab[i]
        res = resTab[i]
        sol = solTab[i]
        (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
        hamiltonian = hamiltonianTab[i]
        hamPhi = hamPhiTab[i]
        hamChi = hamChiTab[i]
        #build difference
        hamiltonian = hamiltonian.-hamiltonian[1]
        hamPhi = hamPhi.-hamPhi[1]
        hamChi = hamChi.-hamChi[1]
        # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian),
#             label=latexstring("|H^{N=",res,"} - H^{0}|"),
#             linewidth=2, linecolor="black", linestyle=:dash
#         )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamPhi),
#             label=latexstring("|H_\\phi^{N=",res,"} - H_\\phi^{0}|"),
            label=latexstring("|H_\\phi - H_\\phi^{0}|"),
            linewidth=2, linecolor=phicolor
        )
        plot!(
            energiesPlot, 
            sol.t, abs.(hamChi),
#             label=latexstring("|H_\\chi^{N=",res,"} - H_\\chi^{0}|"),
            label=latexstring("|H_\\chi - H_\\chi^{0}|"),
            linewidth=2, linecolor=chicolor
        )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian-hamPhi-hamChi),
#             label=latexstring("|V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}|"),
#             linewidth=2, linecolor="orange", linestyle=:dot
#         )
    end

#     for i = 1:(length(resTab)-1)
#         # unpack parameters
#         p = pTab[i]
#         res = resTab[i]
#         sol = solTab[i]
#         (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
#         hamiltonian = hamiltonianTab[i]
#         hamPhi = hamPhiTab[i]
#         hamChi = hamChiTab[i]
#         #build difference
#         hamiltonian = hamiltonian.-hamiltonian[1]
#         hamPhi = hamPhi.-hamPhi[1]
#         hamChi = hamChi.-hamChi[1]
#         # append to plot
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian),
#             label=false,#latexstring("|H^{N=",res,"} - H^{0}|"),
#             linewidth=0.1/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dash
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamPhi),
#             label=false,#latexstring("|H_\\phi^{N=",res,"} - H_\\phi^{0}|"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamChi),
#             label=false,#latexstring("|H_\\chi^{N=",res,"} - H_\\chi^{0}|"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000"
#         )
#         plot!(
#             energiesPlot, 
#             sol.t, abs.(hamiltonian-hamPhi-hamChi),
#             label=false,#latexstring("|V_\\mathrm{int}^{N=",res,"} - V_\\mathrm{int}^{0}|"),
#             linewidth=0.2/(length(resTab) - i), linecolor=colorant"#000000", linestyle=:dot
#         )
#     end

    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "paper_energies_difference_log.png"));
    
    energiesPlot = plot!(
        energiesPlot,
        yscale=:log, 
        ylims=ylimSetAuto,
        yticks=10. .^ range(-11, 11),
        xscale=:log,
        xlim=tlim
    )
    
    # Save the plot to a file (e.g., PNG format)
    savefig(energiesPlot, joinpath(pwd(), dir_path, "paper_energies_difference_loglog.png"));
end



function save_density_plots_paper(sol, p, pInit, dir_path; loss_of_convergence_time=Inf)
    # determine the desired end time
    if !(loss_of_convergence_time isa Number)
        loss_of_convergence_time = sol.t[end]
    end
    loss_of_convergence_time = min(loss_of_convergence_time, sol.t[end])
    
    # unpack parameters
    (dx, mphi2, mchi2, c4, c, epsDiss, Nx, NboundaryPadding) = p
    (a0phi, a0chi, k0phi, k0chi, x0phi, x0chi, offsetphi, offsetchi, aStochastic, mink, maxk, stableRandomSeed) = pInit
    
    # define further parameters
    N = Nx + 2*NboundaryPadding;
    dt = sol.t[2] - sol.t[1];
    Nt = Int(min(length(sol.t), div(loss_of_convergence_time,dt)));

    time = range(0, stop=Nt*dt, length=Nt);
    space = range(0, stop=N*dx, length=N);
    phiField = [sol[:, i][j] for i = 1:Nt, j=1:N]; 
    chiField = [sol[:, i][j] for i = 1:Nt, j=2N+1:3N];
#     phiField = [sol[:, i][j] - offsetphi for i = 1:Nt, j=1:N]; 
#     chiField = [sol[:, i][j] - offsetchi for i = 1:Nt, j=2N+1:3N];
    
    # set density ledgend range
    rangePhi = 2.2; #maximum((maximum(abs.(phiField)),2.1))/2;
    rangeChi = 2.2; #maximum((maximum(abs.(chiField)),2.1))/2;
#     rangePhi = minimum((maximum(abs.(phiField[1:div(Nt, 2),:])),3.5));
#     rangeChi = minimum((maximum(abs.(chiField[1:div(Nt, 2),:])),3.5));
    
    spacingPhi = 10^(floor(log10(rangePhi)))
    rangePhiFloor = floor(rangePhi/spacingPhi)*spacingPhi
    spacingChi = 10^(floor(log10(rangeChi)))
    rangeChiFloor = floor(rangeChi/spacingChi)*spacingChi
    
    # define colors for gradients
    black = RGBA(0, 0, 0, 1)
    white = RGBA(1, 1, 1, 1)
    target_color_phi = RGBA(phicolor.r, phicolor.g, phicolor.b, 1)
    mid_color_phi_1 = RGBA(0.5 * phicolor.r, 0.5 * phicolor.g, 0.5 * phicolor.b, 1)
    mid_color_phi_2 = RGBA(0.5 + 0.5 * phicolor.r, 0.5 + 0.5 * phicolor.g, 0.5 + 0.5 * phicolor.b, 1)
    target_color_chi = RGBA(chicolor.r, chicolor.g, chicolor.b, 1)
    mid_color_chi_1 = RGBA(0.5 * chicolor.r, 0.5 * chicolor.g, 0.5 * chicolor.b, 1)
    mid_color_chi_2 = RGBA(0.5 + 0.5 * chicolor.r, 0.5 + 0.5 * chicolor.g, 0.5 + 0.5 * chicolor.b, 1)
    
    # define color gradients with alpha transparency
    cgradphiAlpha = cgrad([
        RGBA(0, 0, 0, 1),
        RGBA(phicolor.r, phicolor.g, phicolor.b, 0), 
        RGBA(phicolor.r, phicolor.g, phicolor.b, 1)
    ], -1:0.1:1)
    cgradchiAlpha = cgrad([
        RGBA(0, 0, 0, 1),
        RGBA(chicolor.r, chicolor.g, chicolor.b, 0), 
        RGBA(chicolor.r, chicolor.g, chicolor.b, 1)
    ], -1:0.1:1)
    
    # define non-transparent color gradients 
    # (workaround to avoid bad legend behaviour in the single field plots)
    # (for some reason julia's Plots.jl legends to not account for the alpha value)
    cgradphi = cgrad([
            black,
            mid_color_phi_1,
            white,
            mid_color_phi_2,
            target_color_phi
        ], -1:0.1:1)
    cgradchi = cgrad([
            black,
            mid_color_chi_1,
            white,
            mid_color_chi_2,
            target_color_chi
        ], -1:0.1:1)
    
    # generate plot with both fields (and without color bars)
    
    densityplot = plot(
        layout=(1), size = (500,400), framestyle=:box, dpi = 300, 
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{fields}\\;\\;\\phi\\;\\;\\textrm{and}\\;\\;\\chi"),
        xlabel=L"$\textrm{space}\;\;x/L$",
        xticks = 0.2:0.2:0.8, 
        ylabel=L"$\textrm{time}\;\;t/L$",
        left_margin = 40px,
        right_margin = 40px,
        bottom_margin = 0px,
        top_margin = 0px
    )
    heatmap!(
        densityplot,
        space, time, phiField,
        color=cgradphiAlpha,
        cbar=false,
        clim=(-rangePhi,rangePhi)
    )
    heatmap!(
        densityplot,
        space, time, chiField,
        color=cgradchiAlpha,
        cbar=false,
        clim=(-rangeChi,rangeChi)
    )

    savefig(densityplot, joinpath(pwd(), dir_path, "paper_evolution_density.png"));
    
    # generate plot with phi only
    
    # reset density ledend range
    #rangePhi = maximum((maximum(abs.(phiField)),0));
    #rangeChi = maximum((maximum(abs.(chiField)),0));
    
    densityplot = plot(
        layout=(1), size = (500,400), framestyle=:box, dpi = 300, 
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{field}\\;\\;\\phi"),
        xlabel=L"$\textrm{space}\;\;x/L$",
        xticks = 0.2:0.2:0.8, 
        ylabel=L"$\textrm{time}\;\;t/L$",
        left_margin = 40px,
        right_margin = 40px,
        bottom_margin = 0px,
        top_margin = 0px
    )
    heatmap!(
        densityplot,
        space, time, phiField,
        color=cgradphi,
        #cbar=false,
        clim=(-rangePhi,rangePhi),
        colorbar_ticks=-rangePhiFloor:spacingPhi:rangePhiFloor
    )
#     contour!(
#         densityplot,
#         space, time, phiField;
#         levels=10,#Int(rangePhiFloor/spacingPhi), 
#         linewidth=0.1, color=:black
#     )

    savefig(densityplot, joinpath(pwd(), dir_path, "paper_evolution_density_phi.png"));
    
    # generate plot with chi only
    
    densityplot = plot(
        layout=(1), size = (500,400), framestyle=:box, dpi = 300, 
        guidefont = "Computer Modern", tickfont = "Computer Modern",
        xguidefontsize = 16, yguidefontsize = 16, 
        legendfontsize = 16, titlefontsize = 16, tickfontsize = 16,
        #title=latexstring("\\textrm{field}\\;\\;\\chi"),
        xlabel=L"$\textrm{space}\;\;x/L$",
        xticks = 0.2:0.2:0.8, 
        ylabel=L"$\textrm{time}\;\;t/L$",
        left_margin = 40px,
        right_margin = 40px,
        bottom_margin = 0px,
        top_margin = 0px
    )
    heatmap!(
        densityplot,
        space, time, chiField,
        color=cgradchi,
        #cbar=false,
        clim=(-rangeChi,rangeChi),
        colorbar_ticks=-rangeChiFloor:spacingChi:rangeChiFloor
    )
#     contour!(
#         densityplot,
#         space, time, chiField;
#         levels=10,#Int(rangeChiFloor/spacingChi), 
#         linewidth=0.1, color=:black
#     )

    savefig(densityplot, joinpath(pwd(), dir_path, "paper_evolution_density_chi.png"));
    
end
