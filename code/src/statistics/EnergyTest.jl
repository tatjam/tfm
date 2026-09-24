# EnergyTest.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# The energy test compares two sample sets by means of the averages of the
# distances between points of each set, and between points of the two sets.
# The idea is described in detail in
#  "Energy statistics: A class of statistics based on distances" by Székely and Rizzo.
#
# Given sample sets X with n samples and Y with m samples (one sample per COLUMN):
# - (optional) Compute the mean and covariance matrix of X ∪ Y and whiten both
#   sets with them
# - Compute the energy metric
#    E = 2/nm ∑ᵢⱼ ‖Xᵢ-Yⱼ‖ - 1/n² ∑‖Xᵢ-Xⱼ‖ - 1/m² ∑‖Yᵢ-Yⱼ‖
#   and the test statistic  T = (nm)/(n+m) E
# - Estimate the null-distribution of T by randomly splitting X ∪ Y into two
#   sets (sizes n and m) and computing T for each random split
# - The p-value is (1 + #{null samples with T_null ≥ T}) / (1 + #null samples)
#
# To do this relatively fast, we exploit the fact that, if we let D
# be the N×N distance matrix (N = n+m) and w a "label" vector with 1/n on
# the entries belonging to X and -1/m on the entries belonging to Y. Expanding
# the quadratic form wᵀDw gives exactly
#
#    E = -wᵀ D w
#
# This can be evaluated very quickly by BLAS, wayy faster than the direct
# formula, as it can be done with less memory bandwidth (even if, at first glance, it's
# an asymptotically slower algorithm to a clever implementation of the direct formula)

using Statistics: mean, cov
using LinearAlgebra: ldiv!, mul!, dot, Symmetric, cholesky
using Random: shuffle, shuffle!, default_rng, MersenneTwister
using Distances: pairwise, Euclidean

"""
    distance_matrix(Z)

Computes the full N×N (plain `Matrix{Float64}`, both triangles filled, zero
diagonal) matrix D(i, j) = ‖Zᵢ - Zⱼ‖ for the points stored as columns of `Z`.

Despite the distance matrix being symmetric, it's better to not wrap it in
Symmetric() as it adds branches and slows down the SIMD instructions.
"""
function distance_matrix(Z::AbstractMatrix)
    D = pairwise(Euclidean(), convert(Matrix{Float64}, Z), dims = 2)
    @inbounds for i in axes(D, 1)
        D[i, i] = 0.0  # guard against round-off from the BLAS-based Euclidean
    end
    return D
end

"""
    center_distance_matrix!(D)

Double-centers D in place: D := D - rowmeans - colmeansᵀ + grandmean.

The utility of this matrix is evident if we Consider the decomposition of D as
    D = Dᶜ + R + C + G
where
    - Dᶜ is the "centered" matrix
    - R contains, in each element of each row, the mean of said row of D, each ρᵢ
    - C contains, in each element of each col, the mean of said col of D, each χⱼ
    - G is a matrix that contains in each entry the grand mean of D (say γ)

The action of D as a quadratic form in w is then written
    wᵀDw = wᵀ(Dᶜ + R + C + G)w =
         = wᵀ(dw + Rw + Cw + Gw)

Now, note that R = Cᵀ by the fact that our matrix is symmetric. Also, let's assume that
    ∑w = 0 (i.e. its entries sum to zero)

then
    - the action of G on w is 0 by (Gw)ᵢ = γ  ∑w = 0
    - the action of R on w is 0 by (Rw)ᵢ = ρᵢ ∑w = 0
    - the action of C on w is (Cw)ᵢ = ∑ⱼ ρⱼ wⱼ (i.e all rows are the same)

Thus remains wᵀDw = wᵀDᶜw + wᵀ(Cw), the last term is zero by
    wᵀ(Cw) = ∑wᵢ (Cw)ᵢ = (Cw)ᵢ ∑ w = 0 (using the fact all rows of Cw are the same)

Thus we have demonstrated the identity wᵀDw = wᵀDᶜw when ∑w = 0
"""
function center_distance_matrix!(D::AbstractMatrix{Float64})
    # D is symmetric => row means == column means
    r = vec(mean(D, dims = 2))      
    g = mean(r)
    D .= D .- r .- r' .+ g
    return D
end

"""
    centered_distance_matrix(Z, precision = Float64)

Distance matrix of the columns of `Z`, double-centered in Float64 and then
converted to `precision` (Float64 or Float32). See center_distance_matrix! for
more info.
"""
function centered_distance_matrix(
    Z::AbstractMatrix,
    ::Type{T} = Float64,
) where {T<:AbstractFloat}
    D = center_distance_matrix!(distance_matrix(Z))
    return T === Float64 ? D : convert(Matrix{T}, D)
end

# ==============================================================================
# CORE: ENERGY METRIC AS A QUADRATIC FORM  E = -wᵀ D w
# ==============================================================================

"""
    label_weights(n, m, T = Float64)

Label vector of length n+m: 1/n for the first n entries (set X), -1/m for the
remaining m entries (set Y). This ensures ∑w = 0.
"""
function label_weights(n::Int, m::Int, ::Type{T} = Float64) where {T<:AbstractFloat}
    w = fill(-one(T) / m, n + m)
    w[1:n] .= one(T) / n
    return w
end

"""
    energy_metric_from_distance(D, n, m)

Energy metric E for the labelling "first n points are X, last m are Y", given the
(possibly centered) distance matrix D.
"""
function energy_metric_from_distance(D::AbstractMatrix{T}, n::Int, m::Int) where {T}
    w = label_weights(n, m, T)
    return -dot(w, D * w)
end

"""
    energy_statistic_from_distance(D, n, m)

Test statistic T = (nm)/(n+m) E for the labelling "first n points are X".
"""
function energy_statistic_from_distance(D::AbstractMatrix, n::Int, m::Int)
    return (n * m) / (n + m) * energy_metric_from_distance(D, n, m)
end

# Fills each column of W with one random relabelling: `hi` on n random rows,
# `lo` on the other m rows. `perm` is a buffer that must be length n+m.
function _fill_random_labels!(W::AbstractMatrix, perm::Vector{Int}, n::Int, hi, lo, rng)
    for j in axes(W, 2)
        shuffle!(rng, perm)
        col = view(W, :, j)
        fill!(col, lo)
        @inbounds for i = 1:n
            col[perm[i]] = hi
        end
    end
    return W
end

"""
    null_statistics_from_distance!(samples, D, n, m; batch = 128, rng = default_rng())

Fills `samples` with the test statistic of `length(samples)` random relabellings,
`batch` of them per matrix product.
"""
function null_statistics_from_distance!(
    samples::AbstractVector,
    D::AbstractMatrix{T},
    n::Int,
    m::Int;
    batch::Int = 128,
    rng = default_rng(),
) where {T}
    N = n + m
    K = length(samples)
    K == 0 && return samples

    scale = T(n) * T(m) / T(N)
    hi, lo = one(T) / n, -one(T) / m

    B = clamp(batch, 1, K)
    perm = collect(1:N)
    W = Matrix{T}(undef, N, B)
    DW = similar(W)

    k = 1
    while k <= K
        b = min(B, K - k + 1)
        Wb = view(W, :, 1:b)
        DWb = view(DW, :, 1:b)

        _fill_random_labels!(Wb, perm, n, hi, lo, rng)
        mul!(DWb, D, Wb)                       # one GEMM = b permutations

        @inbounds for j = 1:b
            samples[k+j-1] = -scale * dot(view(Wb, :, j), view(DWb, :, j))
        end
        k += b
    end
    return samples
end

"""
    withen_union!(X, Y)

Computes the mean of X ∪ Y and its covariance matrix, and then whitens the points in
both of the sets (in-place) with the whitening associated to said mean and covariance.
The points are given as d×n and d×m matrices (number of points in cols).
"""
function withen_union!(X::AbstractMatrix, Y::AbstractMatrix)
    n, m = size(X, 2), size(Y, 2)
    c = n + m
    cc = c - 1  # Bessel correction

    meanX = mean(X, dims = 2)
    meanY = mean(Y, dims = 2)
    covX = cov(X, dims = 2)
    covY = cov(Y, dims = 2)

    joint_mean = (n * meanX + m * meanY) / c
    meandiff = meanX - meanY
    joint_cov =
        (n - 1) / cc * covX +
        (m - 1) / cc * covY +
        (n * m) / (c * cc) * (meandiff * meandiff')

    L = cholesky(Symmetric(joint_cov)).L

    X .-= joint_mean
    Y .-= joint_mean
    ldiv!(L, X)
    ldiv!(L, Y)

    return X, Y
end

"""
    withen_union(X, Y)

Same as withen_union! but operating on copies of X and Y, returning the whitened samples.
"""
withen_union(X::AbstractMatrix, Y::AbstractMatrix) = withen_union!(copy(X), copy(Y))

"""
    energy_metric(X, Y; precision = Float64)

(Efficiently) computes 2/nm ∑ᵢⱼ ‖Xᵢ-Yⱼ‖ - 1/n² ∑‖Xᵢ-Xⱼ‖ - 1/m²
∑‖Yᵢ-Yⱼ‖ for the points given as d×n and d×m matrices.
(number of points = number of columns).
"""
function energy_metric(X::AbstractMatrix, Y::AbstractMatrix; precision::Type = Float64)
    D = centered_distance_matrix(hcat(X, Y), precision)
    return energy_metric_from_distance(D, size(X, 2), size(Y, 2))
end

"""
    energy_test_statistic(X, Y; precision = Float64)

(Efficiently) computes the test statistic T = (nm)/(n+m) E. X and Y are NOT whitened.
"""
function energy_test_statistic(X::AbstractMatrix, Y::AbstractMatrix; precision::Type = Float64)
    D = centered_distance_matrix(hcat(X, Y), precision)
    return energy_statistic_from_distance(D, size(X, 2), size(Y, 2))
end

"""
    energy_test_null_distribution!(samples, X, Y; batch = 128, rng = default_rng(),
                                   precision = Float64)

(Efficiently) fills `samples` with the test statistic T (same scaling as
`energy_test_statistic`) of random relabellings of X ∪ Y. X and Y are NOT
whitened.
"""
function energy_test_null_distribution!(
    samples::AbstractVector,
    X::AbstractMatrix,
    Y::AbstractMatrix;
    batch::Int = 128,
    rng = default_rng(),
    precision::Type = Float64,
)
    D = centered_distance_matrix(hcat(X, Y), precision)
    return null_statistics_from_distance!(
        samples, D, size(X, 2), size(Y, 2); batch = batch, rng = rng,
    )
end

"""
    energy_test(X, Y, num_samples; whiten = false, batch = 128, rng = default_rng(),
                precision = Float64)

(Efficiently) Performs the energy test on X and Y (one sample per column), using
`num_samples` random permutations to estimate the null-distribution, and returns
the p-value (1 + #{T_null ≥ T}) / (1 + num_samples).

With `whiten = true` the union is whitened first (on copies; X and Y are
untouched). The distance matrix is built once and shared by the observed
statistic and the whole null distribution.
"""
function energy_test(
    X::AbstractMatrix,
    Y::AbstractMatrix,
    num_samples::Int;
    whiten::Bool = false,
    batch::Int = 128,
    rng = default_rng(),
    precision::Type = Float64,
)
    if whiten
        X, Y = withen_union(X, Y)
    end
    n, m = size(X, 2), size(Y, 2)

    D = centered_distance_matrix(hcat(X, Y), precision)
    t = energy_statistic_from_distance(D, n, m)

    null_dist = Vector{Float64}(undef, num_samples)
    null_statistics_from_distance!(null_dist, D, n, m; batch = batch, rng = rng)

    return (1 + count(≥(t), null_dist)) / (1 + num_samples)
end

"""
    sum_pairwise_submatrix(D, idx)

Obtains ∑ D(i, j) where i, j ∈ idx (sorted), exploiting that D is symmetric with
zero diagonal:  ∑ D(i, j) = 2 ∑ D(i, k) for i ∈ idx, k ∈ idx, k < i.

This is a naive algorithm and is fairly slow due to large memory bandwidth needed.
"""
function sum_pairwise_submatrix(D::AbstractMatrix{T}, idx::AbstractVector{Int}) where {T}
    s = zero(T)
    @inbounds for b = 2:length(idx)
        ib = idx[b]
        for a = 1:(b-1)
            s += D[idx[a], ib]
        end
    end
    return 2 * s
end

"""
    energy_metric_on_distance_matrix(D, I, J, D_total)

Energy metric given the RAW (uncentered) distance matrix, where I and J are the
(sorted!) index sets of the two groups and D_total is the sum of all entries of D.

This is a naive algorithm and is fairly slow due to large memory bandwidth needed.
"""
function energy_metric_on_distance_matrix(
    D::AbstractMatrix,
    I::AbstractVector{Int},
    J::AbstractVector{Int},
    D_total,
)
    n, m = length(I), length(J)

    inx = sum_pairwise_submatrix(D, I)
    iny = sum_pairwise_submatrix(D, J)
    cross = 0.5 * (D_total - inx - iny)

    return 2 * cross / (n * m) - inx / (n * n) - iny / (m * m)
end

"""
    energy_metric_reference(X, Y)

Naive evaluation of the energy metric, this is a slow method due to memory bandwidth
"""
function energy_metric_reference(X::AbstractMatrix, Y::AbstractMatrix)
    n, m = size(X, 2), size(Y, 2)
    D = distance_matrix(hcat(X, Y))
    return energy_metric_on_distance_matrix(D, collect(1:n), collect((n+1):(n+m)), sum(D))
end

