%% Nonlinear Analysis of Tower-TMD-FI 
clc; clear;

% Tower
mt = 3.8507e+05;
kt = 2.4671e+06;
zt = 0.01;
wt = sqrt(kt/mt);
Fw = 1e5;
system.mt = mt;

% Damper
mu   = 0.01;
md   = mu*mt;
beta = 0.02;
wr   = 1.706;
zd   = 0.0158;
% zd   = 0;
phi  = 0.5;
b_iner    = beta*mt;
system.beta = beta;
system.phi  = phi;
system.mu   = mu;
system.beta = beta;
system.wt   = wt;
system.zt   = zt;

% [wr, zd] = OptTMDIParams(mu, beta, phi);
% wr = wr*1.1;
system.D_c  = 1.5;
system.D_h  = .8787;
system.r4   = 2.3787;

system.mu_f = 0.0009;
system.rho  = 997;

system.M = [1 + mu + beta*(1-phi)^2    mu + beta*(1-phi);
            mu + beta*(1-phi)              mu + beta];

system.K = [wt^2    0;
            0     mu*wt^2*wr^2];

system.D = [2*wt*zt    0;
            0,   2*mu*wt*wr*zd];

system.n = size(system.M,1);        
n  = size(system.M,1);

%% Compute frequency response using harmonic balance
analysis = 'FRF';

% Analysis parameters
H     = 3;        % harmonic order
N     = 2^9;       % number of time samples per period
delt  = .1;
Om_s = 2;      % start frequency
Om_e = 3;    % end frequency

nu       = Om_s*(1:1:H)';

D_c  = system.D_c;
D_h  = system.D_h;
beta = system.beta;
mu_f = system.mu_f;
rho  = system.rho;
r4   = system.r4;
v = 0;
Fd_norm    = 0.06852*D_c^2*beta/(D_h^(5/2)*sqrt(r4))*v^2 + 35.08*mu_f*beta/(D_h^2*rho)*v;
% length = (b_iner*pi/4*D_h^2)/(rho*(pi/4*D_c^2)^2)
% return

ForceInd = find(ismember(nu, Om_s));
GFInd    = [(1:n)'; kron((ForceInd-1)*n*2+3, ones(2*n,1)) + kron(ones(length(ForceInd),1), (0:2*n-1)')];

system.Fex = zeros((2*H+1)*n,1);
system.Fex(GFInd) = [0, 0, Fw/mt + (phi-1)*Fd_norm, -Fd_norm, 0, 0]';


% Initial guess of solution

Z  = 2*H+1;

Ind = 2:2:Z-1;
vec = ones(Z,1);
vec(Ind) = 0;
vec = [0; kron(nu, [1 1]')].* vec;

A  = kron(eye(Z), system.K) - kron(diag([0; kron(nu, [1; 1])].^2), system.M) + kron(diag(vec(2:end), 1), system.D) - kron(diag(vec(2:end), -1), system.D);

y0 = A\system.Fex;
ds = .05;
Sopt = struct('Dscale',1e-0*ones(length(y0)+1,1),'jac','none');

%% Solve and Plot
% Solve and continue w.r.t. Om
% [Om_HB, a_rms_HB, a_max_HB] = solve_HB(y0,H,system,analysis,Om_s,Om_e,delt,N,ds,Sopt);
% 
% % Plot
% figure(1); hold on;
% plot(Om_HB,a_rms_HB,'linewidth',1); xlim([2 3]); box on;
% xlabel('excitation frequency (\Omega)'); ylabel('response amplitude');
% 
% trapz(Om_HB,a_rms_HB)
% 
% return
%% Optimization

options = optimoptions('fmincon','Algorithm','sqp','OptimalityTolerance', 1e-3,...
            'MaxIterations', 50, 'SpecifyObjectiveGradient', false,...
            'StepTolerance', 1e-3, 'Display', 'iter');

fun = @(x) cost_fun (x,y0,H,system,analysis,Om_s,Om_e,delt,N,ds,Sopt);

A   = [];
b   = [];
Aeq = [];
beq = [];
lb  = [1, .5, .1, -Inf];
ub  = [10, 1.7, 1, Inf];

nonlcontrs = @(x) nonlinear_contraints(x, b_iner, rho);

% lhs_mat = lhsdesign(5,4);
% x0_mat  = [lhs_mat(:,1)*wr, lhs_mat(:,2)*D_c, lhs_mat(:,3)*D_h, lhs_mat(:,4)*r4];
% 
% for i = 1:1
%     x0 = x0_mat(i,:);
%     [x(i,:), fval(i)]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcontrs,options);
% end
x0 = [3 1.5 .3 2];
rng default
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'lb',lb,'ub',ub,'nonlcon',nonlcontrs,'options',options);
x = run(gs,problem)

%% Functions
function [Om_HB, a_rms_HB, a_max_HB] = solve_HB(y0,H,system,analysis,Om_s,Om_e,delt,N,ds,Sopt)


[X_HB,~] = solve_continue(y0,...
    @(X,APDFT) HB_residual(X,APDFT,H,system,analysis),...
    Om_s,Om_e,ds,H,delt,N,Sopt);

% Interpret solver output
Om_HB = X_HB(end,:);
Q_HB = X_HB(1:end-1,:);

% RMS displacement of the 1st degree of freedom (blade)

a_rms_HB = sqrt(sum(Q_HB(3:2:end,:).^2))/sqrt(2);               % this equation is correct as long as mean of displacement(force) is zero

a_max_HB = max(a_rms_HB);
end

function a_max_HB = cost_fun(x,y0,H,system,analysis,Om_s,Om_e,delt,N,ds,Sopt)

mu   = system.mu;
beta = system.beta;
wt   = system.wt;
phi  = system.phi;
zt   = system.zt;
% Optimization variables
wr          = x(1);
system.D_c  = x(2);
system.D_h  = x(3);
system.r4   = x(4);
zd = 0.0001;
% Recalculate system parameters
system.M = [1 + mu + beta*(1-phi)^2    mu + beta*(1-phi);
            mu + beta*(1-phi)              mu + beta];

system.K = [wt^2    0;
            0     mu*wt^2*wr^2];

system.D = [2*wt*zt    0;
            0,   2*mu*wt*wr*zd];

        
[X_HB,~] = solve_continue(y0,...
    @(X,APDFT) HB_residual(X,APDFT,H,system,analysis),...
    Om_s,Om_e,ds,H,delt,N,Sopt);

% Interpret solver output
Q_HB = X_HB(1:end-1,:);

% RMS displacement of the 1st degree of freedom (blade)

a_rms_HB = sqrt(sum(Q_HB(3:2:end,:).^2))/sqrt(2);               % this equation is correct as long as mean of displacement(force) is zero

a_max_HB = max(a_rms_HB);
end

function [c, ceq] = nonlinear_contraints(x, b_iner, rho)
c(1) = (x(2)+x(3))/2 - x(4);
c(2) = x(4) - (x(2)+x(3));
c(3) = (b_iner*pi/4*x(3)^2)/(rho*(pi/4*x(2)^2)^2) - 2;


ceq = [];
end