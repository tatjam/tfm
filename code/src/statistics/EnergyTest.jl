# EnergyTest.jl (c) tatjam 2026
# SPDX-License-Identifier: GPL-3.0-or-later
# ---------------------------------------------
# The energy test is useful to compare the samples of a distribution against a
# set of reference points by means of comparing the averages of the distances
# between points of each sample set, and distances between points between the
# two. The idea is described in detail in
#  "Energy statistics: A class of statistics based on distances" by Székely and Rizzo.
# 
# The idea is, given sample sets X with n samples and Y with m samples we wish
# to compare
# - Compute the mean and covariance matrix of X ∪ Y
# - Withen the samples with this mean and covariance
# - Compute the test statistic
#   T = (nm)/(n+m) (2/nm ∑ᵢⱼ ‖Xᵢ-Yⱼ‖ - 1/n² ∑‖Xᵢ-Xⱼ‖ - 1/m² ∑‖Yᵢ-Yⱼ‖)
# - Estimate the null-distribution of T by randomly sampling X ∪ Y into two
#   classified sets and computing the previous energy metric for the random
#   sets.
# - The p-value of the test is the fraction of the previous samples with
#   T ≥ (the test statistic)
#
# We also implement the distance computation function stand-alone as it's useful
# to see how the methods become incorrect wrt. Monte Carlo as we propagate.

using Statistics: mean, cov
using LinearAlgebra: norm, ldiv!
using Random: shuffle!
using Distances: pairwise, Euclidean

"""
    sum_pairwise_submatrix(D, idx)

Obtains ∑ D(i, j) where i, j ∈ I, which must be sorted, exploiting the fact that D is
Symmetric, and its diagonal is zero, we can compute

∑ D(i, j) = 2∑ D(i, k) for i ∈ I, k ∈ {k ∈ I | k < i}:

  [xxxxx]       [     ]
  [xxxxx]       [x    ]
∑ [xxxxx] = 2 ∑ [xx   ]
  [xxxxx]       [xxx  ]
  [xxxxx]       [xxxx ]
    
"""
function sum_pairwise_submatrix(D, idx)
    len = length(idx)
    s = 0.0

    # By default matrices are column-major, so the innermost loop must iterate
    # along the columns to reduce cache swaps
    @inbounds for b = 2:len
        ib = idx[b]
        @inbounds for a = 1:(b-1)
            ia = idx[a]
            s += D[ia, ib]
        end
    end

    return 2 * s
end

"""
    energy_metric_on_distance_matrix(D, I, J, D_total)
    
Given the distance matrix, where D(i, j) gives distance between points i, j, computes
the energy metric knowing that I is the set of indices which belong to one set and J is
the set of indices which belongs to the other. D_total is the sum of all entries in D.

NOTE: I and J must be sorted before calling this function!
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

    cross /= n * m
    inx /= n * n
    iny /= m * m

    return 2 * cross - inx - iny
end

"""
    distance_matrix(Z)

Computes D(i, j) = ‖Zᵢ - Zⱼ‖ and returns it
    
"""
function distance_matrix(Z::AbstractMatrix)
    return Symmetric(pairwise(Euclidean(), Z, dims = 2))
end

"""
   energy_metric(X, Y) 

Computes 2/nm ∑ᵢⱼ ‖Xᵢ-Yⱼ‖ - 1/n² ∑‖Xᵢ-Xⱼ‖ - 1/m² ∑‖Yᵢ-Yⱼ‖ for the given set of points,
which are given as d×n and d×m matrices (i.e. number of points is the number of columns)
"""
function energy_metric(X::AbstractMatrix, Y::AbstractMatrix)
    n, m = size(X, 2), size(Y, 2)
    Z = hcat(X, Y)
    I = collect(1:n)
    J = collect((n+1):(n+m))
    D = distance_matrix(Z)

    # I and J are sorted by construction
    return energy_metric_on_distance_matrix(D, I, J, sum(D))
end

"""
    withen_union!(X, Y)

Computes the mean of X ∪ Y and its covariance matrix, and then whitens the points in
both of the sets (in-place) with thw whitening associated to said mean and covariance.
Once again, the points are given as d×n and d×m matrices (number of points in cols).
"""
function withen_union!(X::AbstractMatrix, Y::AbstractMatrix)
    n, m = size(X, 2), size(Y, 2)
    c = n + m

    # Bessel correction
    cc = c - 1

    meanX = mean(X, dims = 2)
    meanY = mean(Y, dims = 2)

    covX = cov(X, dims = 2)
    covY = cov(Y, dims = 2)

    mean = (n * meanX + m * meanY) / c
    meandiff = meanX - meanY
    cov =
        (n - 1) / cc * covX +
        (m - 1) / cc * covY +
        (n * m) / (c * cc) * meandiff * meandiff'

    k = cholesky(Symmetric(cov)).L
    ldiv!(X, k, (X .- mean))
    ldiv!(Y, k, (Y .- mean))

    return X, Y
end

"""
    withen_union(X, Y)

Same as whiten_union! but operating on copies of X and Y, returning the whitened samples.  
"""
function withen_union(X::AbstractMatrix, Y::AbstractMatrix)
    return withen_union!(copy(X), copy(Y))
end


"""
   energy_test_statistic(X, Y)

Computes the test statistic
    T = (nm)/(n+m) (2/nm ∑ᵢⱼ ‖Xᵢ-Yⱼ‖ - 1/n² ∑‖Xᵢ-Xⱼ‖ - 1/m² ∑‖Yᵢ-Yⱼ‖)

Note that X and Y are not whitened prior to computing the statistic!
"""
function energy_test_statistic(X::AbstractMatrix, Y::AbstractMatrix)
    n, m = size(X, 2), size(Y, 2)

    return (n * m) / (n + m) * energy_metric(X, Y)
end

"""
    energy_test_null_distribution!(samples, X, Y)

Performs a random sampling of X and Y and computes the energy statistic for each entry in
the samples array, replacing the entry with the resulting statistic.

Once again, X and Y are not whitened before to performing the test.

The random sampling is simply done by stacking both X and Y horizontally (remember,
columns are samples), generating a column array, and repeatedly shuffling it to generate
our partitions. The partitions are just the first n (cols of X) cols, through the indexing
array to the first, and the remainder m (cols of Y) to the second.
"""
function energy_test_null_distribution!(
    samples::AbstractVector,
    X::AbstractMatrix,
    Y::AbstractMatrix,
)
    n, m = size(X, 2), size(Y, 2)
    Z = hcat(X, Y)
    D = distance_matrix(Z)
    D_total = sum(D)

    # {1, 2, ..., n, n+1, n+2, ..., n+m}
    # {X, X, ..., X,   Y,   Y, ...,   Y}
    index_map = collect(1:(n+m))
    for k in eachindex(samples)
        shuffle!(index_map)
        @views I = index_map[1:n]
        @views J = index_map[(n+1):end]
        samples[k] = energy_metric_on_distance_matrix(D, sort(I), sort(J), D_total)
        # @info k
    end
end

"""
    energy_test(X, Y, num_samples)

Performs the energy test to X and Y, using the given number of samples to estimate the
null-distribution, and returning the p-value. The X and Y matrices are NOT whitened, you
must use whiten_union! first if you want to do the whitened energy test.

Each sample is understood to be a column of X or Y.
"""
function energy_test(X::AbstractMatrix, Y::AbstractMatrix, num_samples::Int)
    t = energy_test_statistic(X, Y)
    null_dist = zeros(num_samples)
    energy_test_null_distribution!(null_dist, X, Y)
    return mean(null_dist .≥ t)
end
