function [ xhat_sigma, Wi ] = UKF_sigmaPoints( x_hat, P_xx, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lx      = size(x_hat,1);
mu      = 0;
lambda  = alpha^2 * (lx + mu) - lx;

c       = lambda + lx;

eta     = sqrt(c);
% eta     = sqrt(lx);

%sqrtP   = eta * chol(P_xx,'lower');
[V,D]   = eig(P_xx);
sqrtP   = eta * V*sqrt(abs(D))*V';

% for pp = 1:2*lx+1
%     if pp<lx+1
%         xhat_sigma(:,pp)  = x_hat - sqrtP(:,lx-pp+1);
%     elseif pp==lx+1
%         xhat_sigma(:,pp)  = x_hat;
%     elseif pp>lx+1
%         xhat_sigma(:,pp)  = x_hat + sqrtP(:,pp-lx-1);
%     end
% end

xhat_sigma = [x_hat x_hat-sqrtP x_hat+sqrtP];

Wi = [lambda/c ;
        1/2/c*ones(2*lx,1)];
 


% x_hat(:,ii) = mean(xhat_ii_pp,2);
% y_hat(:,ii) = mean(yhat_ii_pp,2);
% 
% P_xx_iim       = (xhat_ii_pp - x_hat(:,ii)) * (xhat_ii_pp - x_hat(:,ii))' + Q*eye(lx);
% P_yy_iim       = (yhat_ii_pp - y_hat(:,ii)) * (yhat_ii_pp - y_hat(:,ii))' + R*eye(ly);
% P_xx_yy_iim    = (xhat_ii_pp - x_hat(:,ii)) * (yhat_ii_pp - y_hat(:,ii))';
% 
% K_ii    = P_xx_yy_iim * inv(P_yy_iim);
% 
% x_hat(:,ii) = x_hat(:,ii) + K_ii * ( y0(:,ii) - y_hat(:,ii));
% 
% P_xx_iim    = P_xx_iim - K_ii*P_yy_iim*K_ii';

end

