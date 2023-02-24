function fnl = Nonlinear_forces(t, disp, vel, acc, system)

fnl  = zeros(size(disp));
% mt   = system.mt;
D_c  = system.D_c;
D_h  = system.D_h;
beta = system.beta;
mu_f = system.mu_f;
phi  = system.phi;
rho  = system.rho;
r4   = system.r4;

for i = 1:size(fnl,1)
    % Read local variables
    xd_t  = vel(i,1);
    xd_d  = vel(i,2);
    
    v     = xd_d + (1-phi)*xd_t;

    Fd_norm    = 0.06852*D_c^2*beta/(D_h^(5/2)*sqrt(r4))*v*abs(v) + 35.08*mu_f*beta/(D_h^2*rho)*v;

%     Fd_norm    = 0;
        
    fnl(i,1) = (1-phi)*Fd_norm;
    fnl(i,2) = Fd_norm;
    
end