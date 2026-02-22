import Plots

"""
    edge_step_index(odet::OdeState) -> Int

Return a safe edge index for `u_store/ud_store/psi_store` access.
"""
function edge_step_index(odet::OdeState)::Int
    return min(max(odet.step, 1), size(odet.u_store, 4))
end

"""
    mode_block_range(intr::ForceFreeStatesInternal, n_mode::Int) -> UnitRange{Int}

Return linear index range for a given toroidal mode block.
"""
function mode_block_range(intr::ForceFreeStatesInternal, n_mode::Int)::UnitRange{Int}
    n_block = n_mode - intr.nlow + 1
    if n_block < 1 || n_block > intr.npert
        error("n_mode=$n_mode is outside [$(intr.nlow), $(intr.nhigh)]")
    end
    i1 = (n_block - 1) * intr.mpert + 1
    i2 = n_block * intr.mpert
    return i1:i2
end

"""
    evaluate_abc(ffit::FourFitVars, psi::Float64, intr::ForceFreeStatesInternal)

Evaluate full A/B/C matrices at `psi`.
"""
function evaluate_abc(ffit::FourFitVars, psi::Float64, intr::ForceFreeStatesInternal)
    n = intr.numpert_total
    amat = Matrix{ComplexF64}(undef, n, n)
    bmat = Matrix{ComplexF64}(undef, n, n)
    cmat = Matrix{ComplexF64}(undef, n, n)
    ffit.amats(vec(amat), psi; hint=ffit._hint)
    ffit.bmats(vec(bmat), psi; hint=ffit._hint)
    ffit.cmats(vec(cmat), psi; hint=ffit._hint)
    return amat, bmat, cmat
end

"""
    reconstruct_mode_sum(coeffs_m, mvals, theta, n_mode, zeta)

Reconstruct real-space value from Fourier coefficients at one (theta,zeta).
Angles are normalized in [0,1], phase uses exp(2πi(mθ-nζ)).
"""
@inline function reconstruct_mode_sum(
    coeffs_m::AbstractVector{ComplexF64},
    mvals::AbstractVector{Int},
    theta::Float64,
    n_mode::Int,
    zeta::Float64
)::ComplexF64
    val = 0.0 + 0.0im
    @inbounds for (k, m) in pairs(mvals)
        val += coeffs_m[k] * exp(2π * im * (m * theta - n_mode * zeta))
    end
    return val
end

"""
    compute_edge_xi_components(equil, ffit, odet, intr; solution_index=1, n_mode=intr.nlow, zeta=0.0)

Compute edge `xi_psi(θ)`, `xi_s(θ)`, `xi_theta(θ)`, `xi_zeta(θ)` from mode-space data.
`xi_theta` follows the divergence-free ODE described in `papers/DivergenceFreeXi.tex`.
"""
function compute_edge_xi_components(
    equil::Equilibrium.PlasmaEquilibrium,
    ffit::FourFitVars,
    odet::OdeState,
    intr::ForceFreeStatesInternal;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0
)
    step_edge = edge_step_index(odet)
    if step_edge < 2
        error("Need at least 2 saved psi points to estimate ∂ψ(J*xi_psi) at edge.")
    end
    if solution_index < 1 || solution_index > intr.numpert_total
        error("solution_index=$solution_index is outside [1,$(intr.numpert_total)]")
    end

    edge_range = mode_block_range(intr, n_mode)
    mvals = collect(intr.mlow:intr.mhigh)
    psi_edge = odet.psi_store[step_edge]
    psi_prev = odet.psi_store[step_edge-1]
    q_edge = odet.q_store[step_edge]
    dpsi = psi_edge - psi_prev
    if abs(dpsi) < 1e-14
        error("Edge and previous psi are identical; cannot compute radial derivative.")
    end

    # Mode-space vectors for chosen solution column.
    xi_psi_edge_m = @view odet.u_store[edge_range, solution_index, 1, step_edge]
    dxi_psi_edge_m = @view odet.ud_store[edge_range, solution_index, 1, step_edge]
    xi_psi_prev_m = @view odet.u_store[edge_range, solution_index, 1, step_edge-1]

    # Xi_s from A,B,C relation at edge.
    amat, bmat, cmat = evaluate_abc(ffit, psi_edge, intr)
    A = @view amat[edge_range, edge_range]
    B = @view bmat[edge_range, edge_range]
    C = @view cmat[edge_range, edge_range]
    xi_s_edge_m = -(A \ (B * dxi_psi_edge_m + C * xi_psi_edge_m))

    # Stored xi_s for comparison/debug.
    xi_s_edge_stored_m = @view odet.ud_store[edge_range, solution_index, 2, step_edge]

    # Real-space reconstruction on equilibrium theta grid (normalized [0,1]).
    theta = copy(equil.rzphi_ys)
    ntheta = length(theta)
    xi_psi_theta = Vector{ComplexF64}(undef, ntheta)
    xi_s_theta = Vector{ComplexF64}(undef, ntheta)
    xi_psi_prev_theta = Vector{ComplexF64}(undef, ntheta)
    J_edge = Vector{Float64}(undef, ntheta)
    J_prev = Vector{Float64}(undef, ntheta)
    hint_edge = (Ref(1), Ref(1))
    hint_prev = (Ref(1), Ref(1))
    for (j, th) in pairs(theta)
        xi_psi_theta[j] = reconstruct_mode_sum(xi_psi_edge_m, mvals, th, n_mode, zeta)
        xi_s_theta[j] = reconstruct_mode_sum(xi_s_edge_m, mvals, th, n_mode, zeta)
        xi_psi_prev_theta[j] = reconstruct_mode_sum(xi_psi_prev_m, mvals, th, n_mode, zeta)
        J_edge[j] = equil.rzphi_jac((psi_edge, th); hint=hint_edge)
        J_prev[j] = equil.rzphi_jac((psi_prev, th); hint=hint_prev)
    end

    # Source term from divergence-free relation in real space:
    # S = -∂ψ(J xi_psi) - (2π i n / chi1) J xi_s
    chi1 = 2π * equil.psio
    dJxipsi_dpsi = (J_edge .* xi_psi_theta .- J_prev .* xi_psi_prev_theta) ./ dpsi
    S = .-dJxipsi_dpsi .- (2π * im * n_mode / chi1) .* (J_edge .* xi_s_theta)

    # Solve f' - a f = S, f = J xi_theta, with periodic condition over theta∈[0,1].
    a = 2π * im * n_mode * q_edge
    I = zeros(ComplexF64, ntheta)
    for j in 2:ntheta
        h = theta[j] - theta[j-1]
        I[j] = I[j-1] + 0.5 * h * (exp(-a * theta[j-1]) * S[j-1] + exp(-a * theta[j]) * S[j])
    end
    I1 = I[end]
    den = 1 - exp(a)
    if abs(den) < 1e-10
        @warn "Periodic constant denominator near zero at edge (possible rational q≈k/n)."
    end
    f0 = exp(a) * I1 / den
    ftheta = exp.(a .* theta) .* (f0 .+ I)
    xi_theta = ftheta ./ J_edge
    xi_zeta = q_edge .* xi_theta .- xi_s_theta ./ chi1

    return (
        step_edge=step_edge,
        psi_edge=psi_edge,
        psi_prev=psi_prev,
        q_edge=q_edge,
        n_mode=n_mode,
        theta=theta,
        J_edge=J_edge,
        xi_psi=xi_psi_theta,
        xi_s=xi_s_theta,
        xi_theta=xi_theta,
        xi_zeta=xi_zeta,
        xi_s_mode=collect(xi_s_edge_m),
        xi_s_mode_stored=collect(xi_s_edge_stored_m),
        max_mode_xis_diff=maximum(abs.(xi_s_edge_m .- xi_s_edge_stored_m))
    )
end

"""
    print_edge_xi_components(equil, ffit, odet, intr; solution_index=1, n_mode=intr.nlow, zeta=0.0)

Print edge values of `xi_psi`, `xi_theta`, `xi_zeta` for each theta grid point.
Returns the full computed named tuple from `compute_edge_xi_components`.
"""
function print_edge_xi_components(
    equil::Equilibrium.PlasmaEquilibrium,
    ffit::FourFitVars,
    odet::OdeState,
    intr::ForceFreeStatesInternal;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0
)
    data = compute_edge_xi_components(
        equil, ffit, odet, intr;
        solution_index=solution_index, n_mode=n_mode, zeta=zeta
    )
    println("EDGE XI COMPONENTS")
    println("  step=$(data.step_edge), psi_edge=$(data.psi_edge), q_edge=$(data.q_edge), n=$(data.n_mode), zeta=$(zeta)")
    println("  max |Xi_s(ABC) - Xi_s(stored)| in mode space = $(data.max_mode_xis_diff)")
    println("  theta_norm, real(xi_psi), imag(xi_psi), real(xi_theta), imag(xi_theta), real(xi_zeta), imag(xi_zeta)")
    for j in eachindex(data.theta)
        @printf("  %.8f, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e\n",
            data.theta[j],
            real(data.xi_psi[j]), imag(data.xi_psi[j]),
            real(data.xi_theta[j]), imag(data.xi_theta[j]),
            real(data.xi_zeta[j]), imag(data.xi_zeta[j]))
    end
    return data
end

"""
    plot_edge_xi_pretty(data; outpath="edge_xi_vs_theta_pretty.png", dpi=220)

Create a publication-style plot from `compute_edge_xi_components` output.
Returns `(fig, outpath)`.
"""
function plot_edge_xi_pretty(
    data::NamedTuple;
    outpath::String="edge_xi_vs_theta_pretty.png",
    dpi::Int=220
)
    θ = data.theta
    θrad = 2π .* θ

    # Magnitudes and phases
    aψ = abs.(data.xi_psi)
    aθ = abs.(data.xi_theta)
    aζ = abs.(data.xi_zeta)
    pψ = angle.(data.xi_psi)
    pθ = angle.(data.xi_theta)
    pζ = angle.(data.xi_zeta)

    # Shared style
    cψ = :royalblue3
    cθ = :darkorange2
    cζ = :seagreen3
    guide_x = "θ [rad]"

    p1 = Plots.plot(
        θrad, aψ; lw=2.8, color=cψ, label="|xi_psi|",
        xlabel=guide_x, ylabel="Magnitude",
        title="Edge Displacement Magnitudes",
        framestyle=:box, gridalpha=0.20, legend=:topright
    )
    Plots.plot!(p1, θrad, aθ; lw=2.8, color=cθ, label="|xi_theta|")
    Plots.plot!(p1, θrad, aζ; lw=2.8, color=cζ, label="|xi_zeta|")

    p2 = Plots.plot(
        θrad, pψ; lw=2.2, color=cψ, label="arg(xi_psi)",
        xlabel=guide_x, ylabel="Phase [rad]",
        title="Edge Phases (wrapped)",
        framestyle=:box, gridalpha=0.20, legend=:bottomleft
    )
    Plots.plot!(p2, θrad, pθ; lw=2.2, color=cθ, label="arg(xi_theta)")
    Plots.plot!(p2, θrad, pζ; lw=2.2, color=cζ, label="arg(xi_zeta)")
    for xline in (0.0, π, 2π)
        Plots.vline!(p2, [xline]; color=:gray60, lw=1, ls=:dash, label="")
    end

    p3 = Plots.plot(
        θrad, real.(data.xi_psi); lw=2.0, color=cψ, label="Re(xi_psi)",
        xlabel=guide_x, ylabel="Real / Imag",
        title="Real/Imag Components",
        framestyle=:box, gridalpha=0.20, legend=:topright
    )
    Plots.plot!(p3, θrad, imag.(data.xi_psi); lw=2.0, color=cψ, ls=:dash, label="Im(xi_psi)")
    Plots.plot!(p3, θrad, real.(data.xi_theta); lw=2.0, color=cθ, label="Re(xi_theta)")
    Plots.plot!(p3, θrad, imag.(data.xi_theta); lw=2.0, color=cθ, ls=:dash, label="Im(xi_theta)")
    Plots.plot!(p3, θrad, real.(data.xi_zeta); lw=2.0, color=cζ, label="Re(xi_zeta)")
    Plots.plot!(p3, θrad, imag.(data.xi_zeta); lw=2.0, color=cζ, ls=:dash, label="Im(xi_zeta)")

    # Complex-plane trajectories can reveal the "interesting pattern" from phase rotation.
    p4 = Plots.plot(
        real.(data.xi_theta), imag.(data.xi_theta);
        lw=2.0, color=cθ, label="xi_theta(θ)",
        xlabel="Re", ylabel="Im",
        title="Complex-Plane Trajectory",
        framestyle=:box, gridalpha=0.20, aspect_ratio=:equal
    )
    Plots.plot!(p4, real.(data.xi_zeta), imag.(data.xi_zeta); lw=2.0, color=cζ, label="xi_zeta(θ)")
    Plots.scatter!(p4, [real(data.xi_theta[1])], [imag(data.xi_theta[1])]; color=cθ, ms=4, label="start θ=0")

    fig = Plots.plot(
        p1, p2, p3, p4;
        layout=(2, 2), size=(1450, 980), dpi=dpi,
        plot_title="Edge xi Components  (ψ=$(round(data.psi_edge, digits=6)), q=$(round(data.q_edge, digits=4)), n=$(data.n_mode))",
        plot_titlefontsize=12
    )

    Plots.savefig(fig, outpath)
    return fig, outpath
end

"""
    plot_edge_xi_pretty(equil, ffit, odet, intr; kwargs...)

Convenience wrapper: compute edge components then plot.
"""
function plot_edge_xi_pretty(
    equil::Equilibrium.PlasmaEquilibrium,
    ffit::FourFitVars,
    odet::OdeState,
    intr::ForceFreeStatesInternal;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0,
    outpath::String="edge_xi_vs_theta_pretty.png",
    dpi::Int=220
)
    data = compute_edge_xi_components(
        equil, ffit, odet, intr;
        solution_index=solution_index, n_mode=n_mode, zeta=zeta
    )
    fig, saved = plot_edge_xi_pretty(data; outpath=outpath, dpi=dpi)
    return (data=data, fig=fig, outpath=saved)
end

@inline function _trapz_periodic(x::AbstractVector{Float64}, y::AbstractVector)
    n = length(x)
    n == length(y) || error("x/y length mismatch in trapezoidal integration")
    acc = zero(eltype(y))
    @inbounds for i in 1:(n-1)
        acc += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return acc
end

@inline function _periodic_central_derivative(theta::AbstractVector{Float64}, f::AbstractVector)
    n = length(theta)
    n == length(f) || error("theta/f length mismatch in periodic derivative")
    n >= 4 || error("Need at least 4 theta points (including periodic endpoint)")
    nper = n - 1
    df = similar(f)
    @inbounds for j in 1:nper
        if j == 1
            jm = nper
            jp = 2
            θm = theta[jm] - 1.0
            θp = theta[jp]
        elseif j == nper
            jm = nper - 1
            jp = 1
            θm = theta[jm]
            θp = theta[jp] + 1.0
        else
            jm = j - 1
            jp = j + 1
            θm = theta[jm]
            θp = theta[jp]
        end
        df[j] = (f[jp] - f[jm]) / (θp - θm)
    end
    df[end] = df[1]
    return df
end

"""
    compute_edge_pressure_drive_term(equil, odet, intr; solution_index=1, n_mode=intr.nlow, zeta=0.0)

Compute
`-1/(2*μ0) * ∫_0^1 dθ ∫_0^1 dζ J(ψ_edge,θ) * P'(ψ_edge) * |ξ^ψ|^2`
for a selected edge eigen-solution.

Notes:
- `P'` here is `profiles.P_deriv`, i.e. `d(μ0 P)/d(psi_norm)`.
- For a single selected toroidal mode `n`, `|ξ^ψ|^2` is independent of `ζ`,
  so `∫ dζ = 1` and only the `θ` integral is evaluated numerically.
"""
function compute_edge_pressure_drive_term(
    equil::Equilibrium.PlasmaEquilibrium,
    odet::OdeState,
    intr::ForceFreeStatesInternal;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0
)
    step_edge = edge_step_index(odet)
    if step_edge < 1
        error("No saved integration step found in odet")
    end
    if solution_index < 1 || solution_index > intr.numpert_total
        error("solution_index=$solution_index is outside [1,$(intr.numpert_total)]")
    end
    edge_range = mode_block_range(intr, n_mode)
    mvals = collect(intr.mlow:intr.mhigh)
    psi_edge = odet.psi_store[step_edge]
    θ = copy(equil.rzphi_ys)
    ξψ = Vector{ComplexF64}(undef, length(θ))
    J = Vector{Float64}(undef, length(θ))
    hint2d = (Ref(1), Ref(1))
    xi_psi_edge_m = @view odet.u_store[edge_range, solution_index, 1, step_edge]
    for (j, th) in pairs(θ)
        ξψ[j] = reconstruct_mode_sum(xi_psi_edge_m, mvals, th, n_mode, zeta)
        J[j] = equil.rzphi_jac((psi_edge, th); hint=hint2d)
    end
    pprime = equil.profiles.P_deriv(psi_edge; hint=Ref(1))
    int_theta = _trapz_periodic(θ, J .* abs2.(ξψ))
    val = -0.5 / Equilibrium.mu0 * pprime * int_theta

    return (
        value=val,
        psi_edge=psi_edge,
        pprime=pprime,
        int_theta_J_absxi2=int_theta,
        theta=θ,
        xi_psi=ξψ,
        J_edge=J
    )
end

"""
    compute_edge_pressure_drive_term(equil, ffit, odet, intr; solution_index=1, n_mode=intr.nlow, zeta=0.0)

Same as above, but uses provided `ffit` (recommended path).
"""
function compute_edge_pressure_drive_term(
    equil::Equilibrium.PlasmaEquilibrium,
    ffit::FourFitVars,
    odet::OdeState,
    intr::ForceFreeStatesInternal;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0
)
    data = compute_edge_xi_components(
        equil, ffit, odet, intr;
        solution_index=solution_index, n_mode=n_mode, zeta=zeta
    )
    θ = data.theta
    J = data.J_edge
    ξψ = data.xi_psi
    pprime = equil.profiles.P_deriv(data.psi_edge; hint=Ref(1))

    int_theta = _trapz_periodic(θ, J .* abs2.(ξψ))
    val = -0.5 / Equilibrium.mu0 * pprime * int_theta

    return (
        value=val,
        psi_edge=data.psi_edge,
        pprime=pprime,
        int_theta_J_absxi2=int_theta,
        xi_data=data
    )
end

"""
    compute_edge_surface_energy_terms(equil, ffit, odet, intr; solution_index=1, n_mode=intr.nlow, zeta=0.0)

Compute edge surface energy terms:
- `es_discontinuous = -1/(2*mu0) ∫ J * (d(mu0*P)/dpsi) * |xi_psi|^2 dtheta`
- `es_full` from `papers/DivergenceFreeXi.tex` Step 9:
  `-1/2 ∫ J * xi_psi* * B^2 * F dtheta`,
  where `F = div(xi_perp) + 2 xi_perp·kappa`.
"""
function compute_edge_surface_energy_terms(
    equil::Equilibrium.PlasmaEquilibrium,
    ffit::FourFitVars,
    odet::OdeState,
    intr::ForceFreeStatesInternal;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0
)
    data = compute_edge_xi_components(
        equil, ffit, odet, intr;
        solution_index=solution_index, n_mode=n_mode, zeta=zeta
    )

    θ = data.theta
    nθ = length(θ)
    ψ = data.psi_edge
    q = data.q_edge
    chi1 = 2π * equil.psio
    μ0 = Equilibrium.mu0
    J = data.J_edge
    pprime_mu0 = equil.profiles.P_deriv(ψ; hint=Ref(1))

    λ = Vector{ComplexF64}(undef, nθ)                         # λ = xi_parallel/B = (xi·B)/B^2
    xi_perp_dot_gradpsi = Vector{ComplexF64}(undef, nθ)
    xi_perp_dot_gradtheta = Vector{ComplexF64}(undef, nθ)
    div_xi_perp = Vector{ComplexF64}(undef, nθ)
    xi_perp_dot_kappa = Vector{ComplexF64}(undef, nθ)
    F = Vector{ComplexF64}(undef, nθ)
    B2 = Vector{Float64}(undef, nθ)
    dB2_dpsi = Vector{Float64}(undef, nθ)
    dB2_dtheta = Vector{Float64}(undef, nθ)

    hint_r2 = (Ref(1), Ref(1))
    hint_off = (Ref(1), Ref(1))
    hint_nu = (Ref(1), Ref(1))
    hint_B = (Ref(1), Ref(1))

    for (j, th) in pairs(θ)
        r2 = equil.rzphi_rsquared((ψ, th); hint=hint_r2)
        fx1 = equil.rzphi_rsquared((ψ, th); deriv=(1, 0), hint=hint_r2)
        fy1 = equil.rzphi_rsquared((ψ, th); deriv=(0, 1), hint=hint_r2)

        off = equil.rzphi_offset((ψ, th); hint=hint_off)
        fx2 = equil.rzphi_offset((ψ, th); deriv=(1, 0), hint=hint_off)
        fy2 = equil.rzphi_offset((ψ, th); deriv=(0, 1), hint=hint_off)

        fx3 = equil.rzphi_nu((ψ, th); deriv=(1, 0), hint=hint_nu)
        fy3 = equil.rzphi_nu((ψ, th); deriv=(0, 1), hint=hint_nu)

        jac = J[j]
        abs(jac) > 1e-14 || error("Jacobian nearly zero at edge theta index $j")

        rfac = sqrt(max(r2, 0.0))
        rfac > 1e-14 || error("rfac nearly zero at edge theta index $j")
        eta = 2π * (th + off)
        r_major = equil.ro + rfac * cos(eta)

        # Contravariant basis vectors ∇ψ, ∇θ, ∇ζ in cylindrical components.
        v11 = fx1 / (2.0 * rfac * jac)
        v12 = fx2 * 2π * rfac / jac
        v13 = fx3 * r_major / jac
        v21 = fy1 / (2.0 * rfac * jac)
        v22 = (1.0 + fy2) * 2π * rfac / jac
        v23 = fy3 * r_major / jac
        v33 = 2π * r_major / jac

        g_contra = Matrix{Float64}(undef, 3, 3)
        g_contra[1, 1] = v11^2 + v12^2 + v13^2
        g_contra[1, 2] = v11 * v21 + v12 * v22 + v13 * v23
        g_contra[1, 3] = v13 * v33
        g_contra[2, 1] = g_contra[1, 2]
        g_contra[2, 2] = v21^2 + v22^2 + v23^2
        g_contra[2, 3] = v23 * v33
        g_contra[3, 1] = g_contra[1, 3]
        g_contra[3, 2] = g_contra[2, 3]
        g_contra[3, 3] = v33^2
        det(g_contra) > 1e-20 || error("Singular metric tensor at edge theta index $j")
        g_cov = inv(g_contra)

        # Keep B^2 definition consistent with Bal/Fourfit basis:
        # B^2 = chi'^2 * |v2 + q*v3|^2 = chi'^2 * (v21^2 + v22^2 + (v23 + q*v33)^2)
        B2[j] = max(chi1^2 * (v21^2 + v22^2 + (v23 + q * v33)^2), Base.eps(Float64))

        # Use equilibrium spline derivatives for d(B^2)/dpsi and d(B^2)/dtheta.
        B = equil.eqfun_B((ψ, th); hint=hint_B)
        dB_dpsi = equil.eqfun_B((ψ, th); deriv=(1, 0), hint=hint_B)
        dB_dtheta = equil.eqfun_B((ψ, th); deriv=(0, 1), hint=hint_B)
        dB2_dpsi[j] = 2.0 * B * dB_dpsi
        dB2_dtheta[j] = 2.0 * B * dB_dtheta

        Bθ = chi1 / jac
        Bζ = q * Bθ
        b_cov_1 = g_cov[1, 2] * Bθ + g_cov[1, 3] * Bζ
        b_cov_2 = g_cov[2, 2] * Bθ + g_cov[2, 3] * Bζ
        b_cov_3 = g_cov[3, 2] * Bθ + g_cov[3, 3] * Bζ

        ξψ = data.xi_psi[j]
        ξθ = data.xi_theta[j]
        ξζ = data.xi_zeta[j]
        xi_dot_B = ξψ * b_cov_1 + ξθ * b_cov_2 + ξζ * b_cov_3
        λ[j] = xi_dot_B / B2[j]

        ξψ_perp = ξψ
        ξθ_perp = ξθ - λ[j] * Bθ
        ξζ_perp = ξζ - λ[j] * Bζ

        xi_perp_dot_gradpsi[j] = g_cov[1, 1] * ξψ_perp + g_cov[1, 2] * ξθ_perp + g_cov[1, 3] * ξζ_perp
        xi_perp_dot_gradtheta[j] = g_cov[2, 1] * ξψ_perp + g_cov[2, 2] * ξθ_perp + g_cov[2, 3] * ξζ_perp

    end

    dλ_dθ = _periodic_central_derivative(θ, λ)
    div_xi_perp .= .-(chi1 ./ J) .* (dλ_dθ .- (2π * im * n_mode * q) .* λ)
    xi_perp_dot_kappa .= ((pprime_mu0 .+ 0.5 .* dB2_dpsi) .* xi_perp_dot_gradpsi .+ 0.5 .* dB2_dtheta .* xi_perp_dot_gradtheta) ./ B2
    F .= div_xi_perp .+ 2.0 .* xi_perp_dot_kappa

    int_theta_J_absxi2 = _trapz_periodic(θ, J .* abs2.(data.xi_psi))
    es_discontinuous = -0.5 / μ0 * pprime_mu0 * int_theta_J_absxi2
    int_theta_full = _trapz_periodic(θ, J .* conj.(data.xi_psi) .* B2 .* F)
    es_full = -0.5 * int_theta_full

    return (
        psi_edge=ψ,
        q_edge=q,
        n_mode=n_mode,
        pprime_mu0=pprime_mu0,
        es_discontinuous=es_discontinuous,
        es_full=es_full,
        int_theta_J_absxi2=int_theta_J_absxi2,
        int_theta_full=int_theta_full,
        theta=θ,
        xi_data=data,
        lambda_over_B=λ,
        div_xi_perp=div_xi_perp,
        xi_perp_dot_kappa=xi_perp_dot_kappa,
        F=F,
        B2=B2,
        dB2_dpsi=dB2_dpsi,
        dB2_dtheta=dB2_dtheta
    )
end

"""
    print_edge_pressure_drive_comparison(equil, ffit, odet, intr, vac_data; solution_index=1, n_mode=intr.nlow, zeta=0.0)

Print the pressure-drive integral term and compare its order with `ep/ev/et`.
"""
function print_edge_pressure_drive_comparison(
    equil::Equilibrium.PlasmaEquilibrium,
    ffit::FourFitVars,
    odet::OdeState,
    intr::ForceFreeStatesInternal,
    vac_data::VacuumData;
    solution_index::Int=1,
    n_mode::Int=intr.nlow,
    zeta::Float64=0.0
)
    result = compute_edge_surface_energy_terms(
        equil, ffit, odet, intr;
        solution_index=solution_index, n_mode=n_mode, zeta=zeta
    )
    es_discontinuous = result.es_discontinuous
    es_full = result.es_full
    ep = (1 <= solution_index <= length(vac_data.ep)) ? vac_data.ep[solution_index] : NaN + NaN*im
    ev = (1 <= solution_index <= length(vac_data.ev)) ? vac_data.ev[solution_index] : NaN + NaN*im
    et = (1 <= solution_index <= length(vac_data.et)) ? vac_data.et[solution_index] : NaN + NaN*im

    println("EDGE SURFACE ENERGY TERMS")
    @printf("  solution_index = %d, n = %d\n", solution_index, n_mode)
    @printf("  psi_edge       = %.12f\n", result.psi_edge)
    @printf("  q_edge         = %.12f\n", result.q_edge)
    @printf("  pprime_mu0     = %.12e   (d(mu0*P)/d(psi_norm))\n", result.pprime_mu0)
    @printf("  es_discontinuous = %.12e\n", es_discontinuous)
    println("    where es_discontinuous = -1/(2*mu0) * ∫ J * pprime_mu0 * |xi_psi|^2 dtheta")
    @printf("  es_full          = %+.12e %+.12ei\n", real(es_full), imag(es_full))
    println("  --- fixed/free-boundary energies (same solution index) ---")
    println("  --- energy comparison (same solution index) ---")
    @printf("  deltaWp ~ ep   = %+.12e %+.12ei   |ep|=%.12e\n", real(ep), imag(ep), abs(ep))
    @printf("  deltaWv ~ ev   = %+.12e %+.12ei   |ev|=%.12e\n", real(ev), imag(ev), abs(ev))
    @printf("  deltaWt ~ et   = %+.12e %+.12ei   |et|=%.12e\n", real(et), imag(et), abs(et))
    println("  --- real-part comparison (requested) ---")
    @printf("  Re(es_discontinuous) / Re(ep) = %.12e\n", real(es_discontinuous) / max(abs(real(ep)), Base.eps(Float64)))
    @printf("  Re(es_discontinuous) / Re(et) = %.12e\n", real(es_discontinuous) / max(abs(real(et)), Base.eps(Float64)))
    @printf("  Re(es_full)          / Re(ep) = %.12e\n", real(es_full) / max(abs(real(ep)), Base.eps(Float64)))
    @printf("  Re(es_full)          / Re(et) = %.12e\n", real(es_full) / max(abs(real(et)), Base.eps(Float64)))
    return (
        surface_terms=result,
        es_discontinuous=es_discontinuous,
        es_full=es_full,
        ep=ep,
        ev=ev,
        et=et
    )
end

print_edge_surface_energy_comparison(args...; kwargs...) =
    print_edge_pressure_drive_comparison(args...; kwargs...)
