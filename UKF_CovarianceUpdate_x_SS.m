function [ x_prior, P_prior] = UKF_CovarianceUpdate_x_SS( X_k, Q, Wi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lx      = size(X_k,1);
Wd      = diag(Wi);

x_prior = X_k*Wi;
X_tilde = X_k - x_prior;

P_prior   = X_tilde * Wd * X_tilde' + Q(1,1)*eye(lx);


end

