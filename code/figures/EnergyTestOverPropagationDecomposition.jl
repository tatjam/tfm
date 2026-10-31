# EnergyTestOverPropagationDecomposition.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------

using HDF5, Statistics, LinearAlgebra, Printf
using CairoMakie
using OrbitalUncertainty
using SatelliteToolbox

const H5PATH = length(ARGS) ≥ 1 ? ARGS[1] : "comparison.h5"
const OUTDIR = length(ARGS) ≥ 2 ? ARGS[2] : "energy_analysis"

# analyse every STRIDE-th saved step
const STRIDE = 1
# fractions of the time range (exact Shapley, expensive)
const SHAPLEY_TIMES = [0.1, 0.4, 0.7, 1.0]
# reference subsample used for the Shapley subsets
const SHAPLEY_NREF = 1500
# fractions of the time range for witness plots
const WITNESS_TIMES = [0.5, 1.0]
# method indices (1=mc 2=ut 3=stm 4=ut_kep 5=stm_kep)
const WITNESS_METHODS = [2, 3]
# U-statistic (≈0 for identical distributions)
const UNBIASED = true
# clamp for log plots
const FLOOR = 1e-5

const COORDS = ["p", "f", "g", "h", "k", "L"]
const LROW = 6

to_mee(X) = reduce(hcat, [collect(euclid_to_mee(c..., GM_EARTH)) for c in eachcol(X)])
wrapL(row, L0) = L0 .+ mod.(row .- L0 .+ π, 2π) .- π
function to_mee_space(ref, Xs)
    R = to_mee(ref)
    L0 = atan(mean(sin, R[LROW, :]), mean(cos, R[LROW, :]))
    std(R[LROW, :]) > 1 &&
        @warn "Reference L spread > 1 rad: unwrapping around the mean is unreliable" std(
            R[LROW, :],
        )
    R[LROW, :] = wrapL(R[LROW, :], L0)

    Xm = Array{Float64}(undef, 6, size(Xs, 2), size(Xs, 3))
    for m in axes(Xs, 3)
        Mm = to_mee(Matrix(Xs[:, :, m]))
        Mm[LROW, :] = wrapL(Mm[LROW, :], L0)
        Xm[:, :, m] = Mm
    end
    return R, Xm
end

energy2(A, B) = energy_metric(A, B; unbiased = UNBIASED)


# mean distance from each column of Z to all columns of S
function meandist(Z, S)
    d, p = size(Z)
    s = size(S, 2)
    out = zeros(Float64, p)
    Threads.@threads for i = 1:p
        acc = 0.0
        @inbounds for j = 1:s
            q = 0.0
            for k = 1:d
                δ = Z[k, i] - S[k, j]
                q += δ * δ
            end
            acc += sqrt(q)
        end
        out[i] = acc / s
    end
    return out
end

function whiten_info(ref)
    μ = mean(ref, dims = 2)
    L = cholesky(Symmetric(cov(ref, dims = 2))).L
    return (; μ, L, w = X -> L \ (X .- μ))
end

function analyze_time(ref, Xs)
    M = size(Xs, 3)
    wi = whiten_info(ref)
    R = wi.w(ref)

    # per-coordinate standardisation (no decorrelation) -> physical marginals
    μr = vec(mean(ref, dims = 2))
    sd = vec(std(ref, dims = 2))
    stdz = X -> (X .- μr) ./ sd
    Rs = stdz(ref)

    full = zeros(M)
    noloc = zeros(M)
    shape = zeros(M)
    marg = zeros(6, M)
    moff = zeros(M)
    moffphys = zeros(6, M)
    cexc = zeros(M)

    for m = 1:M
        X = Matrix(Xs[:, :, m])
        Xw = wi.w(X)

        full[m] = energy2(Xw, R)

        # 1) remove location error
        mw = vec(mean(Xw, dims = 2))
        noloc[m] = energy2(Xw .- mw, R)

        # 2) also match covariance -> what is left is shape (non-Gaussianity etc.)
        Xsh = whiten_info(Xw).w(Xw)
        shape[m] = energy2(Xsh, R)

        # moments in reference-whitened space
        moff[m] = norm(mw)
        moffphys[:, m] = (wi.L * mw) ./ sd            # mean error per coordinate, in reference σ
        ev = eigen(Symmetric(cov(Xw, dims = 2) - I)).values
        cexc[m] = ev[argmax(abs.(ev))]                # signed largest dispersion mismatch

        # physical-coordinate marginals
        Xstd = stdz(X)
        for j = 1:6
            marg[j, m] = energy2(Xstd[j:j, :], Rs[j:j, :])
        end
    end
    return (; full, noloc, shape, marg, moff, moffphys, cexc)
end

function shapley(ref, Xs; nref = SHAPLEY_NREF)
    ref, Xs = to_mee_space(ref, Xs)
    M = size(Xs, 3)
    n = 6
    μr = vec(mean(ref, dims = 2))
    sd = vec(std(ref, dims = 2))
    stdz = X -> (X .- μr) ./ sd
    Rs = stdz(ref[:, 1:min(nref, size(ref, 2))])
    Xstd = [stdz(Matrix(Xs[:, :, m])) for m = 1:M]

    vals = zeros(M, 2^n)
    for mask = 1:(2^n-1)
        idx = [j for j = 1:n if mask & (1 << (j - 1)) != 0]
        Rm = Rs[idx, :]
        for m = 1:M
            vals[m, mask+1] = energy2(Xstd[m][idx, :], Rm)
        end
    end

    φ = zeros(M, n)
    for j = 1:n, mask = 0:(2^n-1)
        mask & (1 << (j - 1)) != 0 && continue
        s = count_ones(mask)
        wgt = factorial(s) * factorial(n - s - 1) / factorial(n)
        for m = 1:M
            φ[m, j] += wgt * (vals[m, (mask|(1<<(j-1)))+1] - vals[m, mask+1])
        end
    end
    return φ        # rows sum to v(full) (efficiency)
end

# f(z) = E|z - ref| - E|z - X|   (>0: z lies where the method has excess mass)
function witness_plot(ref, X, title, path)
    N = size(ref, 2)
    μr = vec(mean(ref, dims = 2))
    sd = vec(std(ref, dims = 2))
    stdz = Z -> (Z .- μr) ./ sd
    Rs, Xst = stdz(ref), stdz(X)

    h = N ÷ 2                       # disjoint halves: evaluation points vs reference
    Rref = Rs[:, 1:min(h, 1500)]
    Zr = Rs[:, (h+1):min(N, h+800)]
    wr = meandist(Zr, Rref) .- meandist(Zr, Xst)
    wx = meandist(Xst, Rref) .- meandist(Xst, Xst)

    # principal axes of the reference (standardised coordinates)
    F = eigen(Symmetric(cov(Rs, dims = 2)))
    V = F.vectors[:, sortperm(F.values, rev = true)][:, 1:3]
    Pr, Px = V' * Zr, V' * Xst

    lim = maximum(abs, vcat(wr, wx))
    fig = Figure(size = (1250, 760))
    Label(fig[0, 1:3], title; fontsize = 18)
    pairs = [(1, 2), (1, 3), (2, 3)]
    local sc
    for (c, (a, b)) in enumerate(pairs)
        ax1 = Axis(
            fig[1, c];
            xlabel = "PC$a [σ]",
            ylabel = "PC$b [σ]",
            title = c == 1 ? "reference points" : "",
        )
        sc = scatter!(
            ax1,
            Pr[a, :],
            Pr[b, :];
            color = wr,
            colormap = :balance,
            colorrange = (-lim, lim),
            markersize = 6,
        )
        ax2 = Axis(
            fig[2, c];
            xlabel = "PC$a [σ]",
            ylabel = "PC$b [σ]",
            title = c == 1 ? "method points" : "",
        )
        scatter!(
            ax2,
            Px[a, :],
            Px[b, :];
            color = wx,
            colormap = :balance,
            colorrange = (-lim, lim),
            markersize = 6,
        )
    end
    Colorbar(fig[1:2, 4], sc; label = "witness  E|z−ref| − E|z−method|")
    save(path, fig)
end

function main()
    mkpath(OUTDIR)
    h5open(H5PATH, "r") do f
        Δt = Float64(attrs(f)["delta_t"])
        names = String.(split(attrs(f)["methods"], ","))
        M = length(names)
        nsteps = size(f["energies"], 2)
        dtruth, dsamp = f["truth"], f["samples"]

        # skip trailing steps that were never written (interrupted run)
        last = nsteps
        while last > 1 && !any(!iszero, dtruth[:, 1:5, last])
            last -= 1
        end
        idxs = unique(vcat(collect(1:STRIDE:last), last))
        T = length(idxs)
        # NOTE: generator writes index i at t = Δt*i (not ts[i]); we follow that.
        td = Δt .* idxs ./ 86400
        @info "Analysing" file = H5PATH steps_written = last analysed = T methods = names

        # sanity check: V-statistic vs the saved Float32 energies (tells you the convention)
        let i = idxs[cld(T, 2)]
            ref = Float64.(dtruth[:, :, i])
            X = Float64.(dsamp[:, :, 1, i])
            w = whiten_info(ref).w
            d2v = energy_metric(w(X), w(ref))
            @printf(
                "Sanity (mc, step %d): saved = %.6g | V-stat D² = %.6g | sqrt = %.6g\n",
                i,
                f["energies"][1, i],
                d2v,
                sqrt(max(d2v, 0))
            )
        end

        full = zeros(M, T)
        noloc = zeros(M, T)
        shape = zeros(M, T)
        marg = zeros(6, M, T)
        moff = zeros(M, T)
        moffphys = zeros(6, M, T)
        cexc = zeros(M, T)

        for (p, i) in enumerate(idxs)
            ref = Float64.(dtruth[:, :, i])
            Xs = Float64.(dsamp[:, :, :, i])
            r = analyze_time(ref, Xs)
            full[:, p] = r.full
            noloc[:, p] = r.noloc
            shape[:, p] = r.shape
            marg[:, :, p] = r.marg
            moff[:, p] = r.moff
            moffphys[:, :, p] = r.moffphys
            cexc[:, p] = r.cexc
            p % 10 == 0 && @info "Decomposition" p T
        end

        # ---- Shapley at selected times
        pos(fr) = clamp(round(Int, fr * T), 1, T)
        spos = unique(pos.(SHAPLEY_TIMES))
        shap = Dict{Int,Matrix{Float64}}()
        for p in spos
            @info "Shapley" t_days = td[p]
            i = idxs[p]
            shap[p] = shapley(Float64.(dtruth[:, :, i]), Float64.(dsamp[:, :, :, i]))
        end

        # ---- CSV
        open(joinpath(OUTDIR, "decomposition.csv"), "w") do io
            println(
                io,
                "t_days,method,full,no_loc,shape,mean_err_norm,cov_excess,",
                join(("marg_" .* COORDS), ","),
                ",",
                join(("meanerr_sigma_" .* COORDS), ","),
            )
            for p = 1:T, m = 1:M
                println(
                    io,
                    join(
                        [
                            td[p];
                            names[m];
                            full[m, p];
                            noloc[m, p];
                            shape[m, p];
                            moff[m, p];
                            cexc[m, p];
                            marg[:, m, p];
                            moffphys[:, m, p]
                        ],
                        ",",
                    ),
                )
            end
        end

        # ---- console summary
        for p in spos
            @printf("\n=== t = %.2f d ===\n", td[p])
            @printf(
                "%-9s %11s %11s %11s %9s %9s\n",
                "method",
                "D²",
                "no-mean",
                "shape-only",
                "|mean|σ",
                "cov-exc"
            )
            for m = 1:M
                @printf(
                    "%-9s %11.4g %11.4g %11.4g %9.3g %9.3g\n",
                    names[m],
                    full[m, p],
                    noloc[m, p],
                    shape[m, p],
                    moff[m, p],
                    cexc[m, p]
                )
            end
            @printf(
                "Shapley (standardised coords, sums to D²-full of the subsample)\n%-9s %s\n",
                "",
                join(lpad.(COORDS, 10))
            )
            for m = 1:M
                @printf(
                    "%-9s %s\n",
                    names[m],
                    join([@sprintf("%10.3g", v) for v in shap[p][m, :]])
                )
            end
        end

        # ---- figures
        cols = Makie.wong_colors()

        fig = Figure(size = (900, 520))
        ax = Axis(
            fig[1, 1];
            xlabel = "t [days]",
            ylabel = "D² (whitened, U-stat)",
            yscale = log10,
            title = "Energy distance to ground truth",
        )
        for m = 1:M
            lines!(ax, td, max.(full[m, :], FLOOR); label = names[m], color = cols[m])
        end
        axislegend(ax; position = :lt)
        save(joinpath(OUTDIR, "energy_overview.png"), fig)

        fig = Figure(size = (1350, 750))
        for m = 1:M
            ax = Axis(
                fig[div(m-1, 3)+1, mod(m-1, 3)+1];
                title = names[m],
                xlabel = "t [days]",
                ylabel = "D²",
                yscale = log10,
            )
            lines!(ax, td, max.(full[m, :], FLOOR); label = "full", color = :black)
            lines!(
                ax,
                td,
                max.(noloc[m, :], FLOOR);
                label = "mean removed",
                color = cols[2],
            )
            lines!(
                ax,
                td,
                max.(shape[m, :], FLOOR);
                label = "mean+cov matched (shape)",
                color = cols[3],
            )
            lines!(
                ax,
                td,
                max.(full[1, :], FLOOR);
                label = "MC floor",
                color = :gray,
                linestyle = :dash,
            )
        end
        Legend(fig[2, 3], first(filter(b -> b isa Axis, fig.content)); framevisible = false)
        save(joinpath(OUTDIR, "decomposition.png"), fig)

        fig = Figure(size = (1350, 750))
        for m = 1:M
            ax = Axis(
                fig[div(m-1, 3)+1, mod(m-1, 3)+1];
                title = names[m] * " – marginals",
                xlabel = "t [days]",
                ylabel = "1D D² (standardised)",
                yscale = log10,
            )
            for j = 1:6
                lines!(
                    ax,
                    td,
                    max.(marg[j, m, :], FLOOR);
                    label = COORDS[j],
                    color = cols[mod1(j, 7)],
                )
            end
        end
        Legend(fig[2, 3], first(filter(b -> b isa Axis, fig.content)); framevisible = false)
        save(joinpath(OUTDIR, "marginals.png"), fig)

        fig = Figure(size = (1000, 720))
        ax1 = Axis(
            fig[1, 1];
            ylabel = "‖mean error‖ (ref-whitened)",
            yscale = log10,
            title = "Location error",
        )
        ax2 = Axis(
            fig[2, 1];
            xlabel = "t [days]",
            ylabel = "largest |eig(cov − I)| (signed)",
            title = "Dispersion mismatch (>0 over-dispersed, <0 under-dispersed)",
            yscale = Makie.Symlog10(1e-2),
        )
        for m = 1:M
            lines!(ax1, td, max.(moff[m, :], FLOOR); label = names[m], color = cols[m])
            lines!(ax2, td, cexc[m, :]; color = cols[m])
        end
        axislegend(ax1; position = :lt)
        save(joinpath(OUTDIR, "moments.png"), fig)

        fig = Figure(size = (400 * length(spos), 450))
        for (c, p) in enumerate(spos)
            ax = Axis(
                fig[1, c];
                title = @sprintf("Shapley, t = %.1f d", td[p]),
                xticks = (1:6, COORDS),
                ylabel = c == 1 ? "contribution to D²" : "",
            )
            w = 0.8 / M
            for m = 1:M
                barplot!(
                    ax,
                    (1:6) .+ (m - (M + 1) / 2) * w,
                    shap[p][m, :];
                    width = w,
                    color = cols[m],
                    label = names[m],
                )
            end
            c == 1 && axislegend(ax; position = :lt)
        end
        save(joinpath(OUTDIR, "shapley.png"), fig)

        # ---- witness plots
        for fr in WITNESS_TIMES, m in WITNESS_METHODS
            p = pos(fr)
            i = idxs[p]
            ref = Float64.(dtruth[:, :, i])
            X = Float64.(dsamp[:, :, m, i])
            witness_plot(
                ref,
                X,
                @sprintf("Witness: %s vs truth, t = %.1f d", names[m], td[p]),
                joinpath(OUTDIR, @sprintf("witness_%s_t%.1fd.png", names[m], td[p])),
            )
        end
    end
    @info "Done" OUTDIR
end

main()
