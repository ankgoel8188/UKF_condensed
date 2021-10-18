clc
clear all
close all
rand('state',1)
randn('state',1)
% LoadFigurePrintingProperties
set(0, 'DefaultLineLineWidth', 1);


n = 3;

A = .3*rand(n,n);
B = rand(n,1);
C = rand(1,n);
Q = 2;1e-9;
R = 1;


a   = 0.8 + 0.6i;

Num = poly([0.3 0.2]);
Den = poly([a a' -0.5]);
% Den = poly([a a' 0.5]);
% G = tf(Num, Den,1)
% pzmap(G)
[A,B,C,D] = tf2ss(Num,Den);
B   = round(10*rand(size(A,2), 1))/10;
C   = round(10*randn(1,size(A,1)))/10;


lx  = size(A,1);
lu  = size(B,2);
lz  = size(C,1);
ly  = lz;

Q   = 1e-8*eye(lx);
R   = 1e-6*eye(ly);
Gamma = eye(n);

System.A = A;
System.B = B;
System.C = C;
System.Q = Q;
System.R = R;

I_lx = eye(lx);

N = 20000;
N = 5000;
N = 50;
%%





rand('state',1)
randn('state',1)



x       = zeros(lx,N);
y       = zeros(ly,N);
u       = zeros(1,N);

xhat_OP     = 0*x;
xhat_OF     = 0*x;
xhat_UKF    = 0*x;

lmu     = size(Gamma,2);
leta    = lmu;
mu      = zeros(lmu,N);
z       = zeros(ly,N);

P_OP   = 1e1*eye(n);
P_OF   = P_OP; 10*eye(n);
P_UKF  = 1e+1*eye(lmu);
P_UKF_SS  = P_OP;
P_UKF_FS  = P_OP;

K_OP    = zeros(lmu, N);
K_OF    = K_OP;
K_UKF   = K_OP;

J_OP    = zeros(1, N);
J_OF    = J_OP;
J_UKF   = J_OP;

err_OP    = zeros(1, N);
err_OF    = err_OP;
err_UKF   = err_OP;

x(:,1)  = [1 1 1]'+sqrt(Q)*randn(lx,1);
y(:,1)  = C*x(:,1)+sqrt(R)*randn(ly,1);


for ii = 1:N
    
    u(:,ii)     = 0*randn;
    x(:,ii+1)   = A*x(:,ii) + B*u(:,ii) + sqrt(Q)*randn(lx,1);
    y(:,ii+1)   = C*x(:,ii+1) + sqrt(R)*randn(ly,1);
    
    % KF
    [K_OF(:,ii),P_OF,J_OF(:,ii)] = Compute_KPJ(System, P_OF, 'OF' );
    xbar_OF         = A*xhat_OF(:,ii) + B*u(:,ii);
    xhat_OF(:,ii+1) = xbar_OF + Gamma*K_OF(:,ii) * ( y(:,ii+1) - C*xbar_OF);
    err_OF(ii)      = norm(xhat_OF(:,ii)-x(:,ii))^2;
    
    
    % UKF
    system = System;
    system.uk = u(:,ii);
    UKF_data.alpha  = 1.2;
    UKF_data.ly     = ly;
    [xhat_UKF(:,ii+1), P_UKF, K_UKF(:,ii) ] = UKF_propagate_est_cov(xhat_UKF(:,ii), P_UKF, y(:,ii+1), @LinSys_dynamics, @LinSys_output, system, UKF_data);
    
    err_UKF(ii) = norm(xhat_UKF(:,ii)-x(:,ii))^2;
    J_UKF(:,ii) = trace(P_UKF);
    
end


figure(3)
row = 1;
green = [0    0.4470    0.7410]; [.75 .75 .75];
subplot(row,2,1)
stairs(err_OF', 'linewidth',2)
hold on; axis tight; grid on
stairs(J_OF, 'linewidth',2)
set(gca,'yscale', 'log')
legend('|e_k|', 'tr P_k')

subplot(row,2,2)
stairs(err_UKF','linewidth',2)
hold on; axis tight; grid on
stairs(J_UKF,'linewidth',2)
set(gca,'yscale', 'log')
legend('|e_k|', 'tr P_k')

