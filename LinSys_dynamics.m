function xkp1 = LinSys_dynamics(xk,system)
   xkp1 = system.A*xk + system.B*system.uk; 
end