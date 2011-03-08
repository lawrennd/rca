function X = mgd_sample(rmean,covariance)
% Generates samples from a Multivariate Gaussian Distribution
%
% X = mgd(mu, covariance) returns N Gaussian-distributed samples in D
% dimensions. mu is a N-by-D matrix and the ith row contains the mean of
% the Gaussian distribution from which the ith row of X is sampled.
% covariance is a D-by-D positive-definite matrix with the diagonal terms
% being the variances of the respective components of x, and the
% off-diagonal terms being the covariances between the different dimensions
% of x. mgd_sample differs from randn in that the mean and covariance can
% be specified.
%
% Example: Generate 50 samples from a 2 dimensional multivariate Gaussian
% distribution, with a mean of [4 5] and a covariance matrix of [9 0;0 9].
% X = mgd(ones(50,1)*[4 5], [9 0;0 9])
%
% mgd_sample is inspired by 'mgd' at
% http://www.mathworks.com/matlabcentral/fx_files/5984/1/mgd.m
% Alfredo Kalaitzis, 2009

% parse input parameters
[N,D]=size(rmean);
[rowsv,colsv]=size(covariance);

if rowsv ~= colsv
    error('Covariance matrix should be square')
end
if D ~= colsv
    error('The dimension of the requested mean is not equal to the dimensions of the covariance')
end

X = randn(N,D);  % generate the samples using built in Matlab function
xmean = mean(X); % calculate the sample mean

% remove any mean from the Matlab generated numbers
% as N increases the Matlab mean approaches zero
X = X - ones(N,1)*xmean;

% make all Gaussian distributions have the covariance specified by the
% user. compute the Cholesky decomposition of the given covariance matrix.
% use a different method depending on wheter or not the covariance is
% positive definite. This is used to scale each component of a sample by
% the square root of the variance.
[R,p] = chol(covariance); % R upper triangular, p flags pos-definiteness
if p > 0
    X = X * sqrtm(covariance);
else
    X = X * R;
end

% add back on the distribution means
X = X + rmean;
end

