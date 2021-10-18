function [P_post,K_k, y_pr] = UKF_CovarianceUpdate_xy_SS2( X_k, Y_k, P_prior, R, Wi)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ly      = size(Y_k,1);
Wd      = diag(Wi);

x_prior = X_k*Wi;
y_pr    = Y_k*Wi;

X_tilde = X_k - x_prior;
Y_tilde = Y_k - y_pr;


P_yy_k  = Y_tilde * Wd * Y_tilde' + R*eye(ly);
P_xy_k  = X_tilde * Wd * Y_tilde';


K_k     = P_xy_k * pinv(P_yy_k);
P_post  = P_prior - K_k * P_xy_k';






% u       = K_k * ( y_k - y_pr);








end

