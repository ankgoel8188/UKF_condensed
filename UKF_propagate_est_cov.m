function [x_post, P_post, K_kp1 ] = UKF_propagate_est_cov(x_k, P_k, y_kp1, dynamics, output, system, UKF_data)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here


% lx = size(X_k,1);

[X_k , Wi]= UKF_sigmaPoints( x_k, P_k, UKF_data.alpha);

X_kp1 = X_k*0;
for ii = 1:size(X_k,2)    
    X_kp1(:,ii) = dynamics(X_k(:,ii), system);
end

[x_prior, P_prior] = UKF_CovarianceUpdate_x_SS( X_kp1, system.Q, Wi );
    

[X_kp1_op , Wi]= UKF_sigmaPoints( x_prior, P_prior, UKF_data.alpha);


Y_kp1 = zeros(UKF_data.ly, size(X_k,2));
for ii = 1:size(X_k,2)    
    Y_kp1(:,ii) = output(X_kp1_op(:,ii), system);
end


[P_post,K_kp1, y_pr] = UKF_CovarianceUpdate_xy_SS2( X_kp1_op, Y_kp1, P_prior, system.R, Wi);


x_post = x_prior + K_kp1 * (y_kp1 - y_pr);


end




