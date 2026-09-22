# GVM.jl (c) tatjam 2026
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

"""
   energy_metric(X, Y) 

Computes 2/nm ∑ᵢⱼ ‖Xᵢ-Yⱼ‖ - 1/n² ∑‖Xᵢ-Xⱼ‖ - 1/m² ∑‖Yᵢ-Yⱼ‖ for the given set of points,
which are given as d×n and d×m matrices (i.e. number of points is the number of columns)
"""
function energy_metric(X::AbstractMatrix, Y::AbstractMatrix)
    n, m = size(X, 2), size(Y, 2)
    cross = mean(norm(X[:,i] - Y[:,j]) for i in 1:n, j in 1:m)
    inx = mean(norm(X[:,i] - Y[:,j]) for i in 1:n, j in 1:n)
    iny = mean(norm(X[:,i] - Y[:,j]) for i in 1:n, j in 1:n)
    return 2 * cross - inx - iny
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

    # For covariances we use the Bessel correction
    cc = c - 2

    meanX = mean(X, dims=2)
    meanY = mean(Y, dims=2)

    covX = cov(X, dims=2)
    covY = cov(Y, dims=2)

    mean = (n * meanX + m * meanY) / c
    meandiff = meanX - meanY
    cov = (n - 1) / cc * covX + (m - 1) / cc * covY + (n * m) / (c * cc) * meandiff * meandiff'

    k = cholesky(Symmetric(cov)).L
    X .= k \ (X .- mean) 
    Y .= k \ (Y .- mean)

    return X, Y
end

"""
    withen_union(X, Y)

Same as whiten_union! but operating on copies of X and Y, returning the whitened sample sets.  
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
    energy_test(X, Y)

Performs the energy test to X and Y, using the given number of samples to estimate the
null-distribution, and returning the p-value. The X and Y matrices are NOT whitened, you
must use whiten_union! first if you want to do the whitened energy test.
"""
function energy_test(X::AbstractMatrix, Y::AbstractMatrix, num_samples::Int)
    t = energy_test_statistic(X, Y)

    
end

