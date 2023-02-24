function [Gamma, t] = APDFT_fun(delt, N, nu)

% Check deltat 
% DT(:,1) = 2*pi*(1:10)'./(kron(nu, ones(length(nu),1)) + kron(ones(length(nu),1), nu));
% DT(:,2) = 2*pi*(1:10)'./(kron(nu, ones(length(nu),1)) - kron(ones(length(nu),1), nu));

M = length(nu);
Gamma = ones(N+1, 2*M+1);
t = 0:delt:N*delt;

for i = 1:N+1
   for j = 2:2:2*M+1

       Gamma(i, j:j+1) = [cos(nu(j/2)*(i-1)*delt) sin(nu(j/2)*(i-1)*delt)];
       
   end
end
