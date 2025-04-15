phicolor = colorant"#c51b8a";
chicolor = colorant"#00FFFF";




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
    phiField = [sol[:, i][j] - offsetphi for i = 1:Nt, j=1:N]; 
    chiField = [sol[:, i][j] - offsetchi for i = 1:Nt, j=2N+1:3N];
    
    # set density ledgend range
    rangePhi = 1.05;
    rangeChi = 1.05;
    
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
        layout=(1), size = (400,600), framestyle=:box, dpi = 400, 
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
        layout=(1), size = (400,600), framestyle=:box, dpi = 400, 
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
    contour!(
        densityplot,
        space, time, phiField;
        levels=11,#Int(rangePhiFloor/spacingPhi), 
        linewidth=0.1, color=:black
    )

    savefig(densityplot, joinpath(pwd(), dir_path, "paper_legend_phi.png"));
    
    # generate plot with chi only
    
    densityplot = plot(
        layout=(1), size = (400,600), framestyle=:box, dpi = 400, 
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
    contour!(
        densityplot,
        space, time, chiField;
        levels=11,#Int(rangeChiFloor/spacingChi), 
        linewidth=0.1, color=:black
    )

    savefig(densityplot, joinpath(pwd(), dir_path, "paper_legend_chi.png"));
    
end
