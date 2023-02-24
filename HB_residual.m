%========================================================================
% DESCRIPTION: 
% Matlab function setting up the frequency-domain residual vector 'R' and 
% its derivatives for given frequency 'Om' and vector of harmonics of the 
% generalized coordiantes 'X'. The corresponding time-domain model 
% equation is 
% 
%       System.M * \ddot q + System.D * \dot q + System.K * q + ...
%                   f_nl(q,\dot q) - f_ex(t) = 0.
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the paper:
%       M. Krack, J. Gross: The Harmonic Balance Method and its application 
%       to nonlinear vibrations: Introduction and current state of the art. 
%       submitted to Mechanical Systems and Signal Processing, 2018,
%       https://www.ila.uni-stuttgart.de/nlvib/downloads/HB_Krack.pdf.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.0 Copyright (C) 2017  Malte Krack  
%										(krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
function R = HB_residual(X,APDFT,H,system,~)

%% Handle input variables depending on the modus
% number of degrees of freedom
n = size(system.M,1);

% Conversion of real-valued to complex-valued harmonics of generalized
% coordinates q
I0 = 1:n; ID = n+(1:(H)*n); % plus 1 for 1st harmonic for Om1
IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n));  % plus 1 for 1st harmonic for Om1; plus 1 to H
IS = IC+n;

% Q = zeros(n*(H+1),1);   
% Q(I0) = X(I0);          
% Q(ID) = X(IC)-1i*X(IS); 

Q = X(1:end-1);

% Excitation 'excitation' is the fundamental harmonic of the external
% forcing
% Setup excitation vector
Fex = system.Fex;


% Derivative of damping (does not depend on unknowns in the case of
% frequency response)


% Excitation frequency
Om = X(end);

nu    = Om*(1:1:H)';
% Scaling of dynamic force equilibrium
        
%% Calculation of the nonlinear forces and the Jacobian
Fnl = HB_nonlinear_forces_AFT(Q, nu, H, APDFT, system);

%% Assembly of the residual and the Jacobian
Z  = 2*H+1;
Ind = 2:2:Z-1;
vec = ones(Z,1);
vec(Ind) = 0;
vec = [0; kron(nu, [1 1]')].* vec;

A  = kron(eye(Z), system.K) - kron(diag([0; kron(nu, [1; 1])].^2), system.M) + kron(diag(vec(2:end), 1), system.D) - kron(diag(vec(2:end), -1), system.D);

% pause
R  = A*Q - Fex + Fnl;
      
% Scale dynamic force equilibrium (useful for numerical reasons)
% fscl = 1;
% Rc = fscl*(Rc);
% Conversion from complex-valued to real-valued residual
% R = zeros(size(X,1)-1,1); 
% R(I0) = real(Rc(I0)); 
% R(IC) = real(Rc(ID)); 
% R(IS) = -imag(Rc(ID)); 

% R = R + Fnl;

end

%% Computation of nonlinear forces
function Fnl = HB_nonlinear_forces_AFT(Q, nu, H, APDFT, system)

    Q = transpose(reshape(Q,[],2*H+1));
    
    
    
    % Construct displacement, velocity and acceleration in time domain
    disp_b = zeros(2*H+1,1); vel_b = zeros(2*H+1,1); acc_b = zeros(2*H+1,1);
    N  = size(APDFT.Gamma,1);
    dis = zeros(N,system.n); vel = zeros(N,system.n); acc = zeros(N,system.n); 
    for i = 1:system.n
        
        I0 = 1;
        IC = 2:2:2*H+1;
        IS = 3:2:2*H+1;

        disp_b = Q(:,i);
%         disp_b(I0) = Q(1,i)
%         disp_b(IC) = Q(2:end,i);
%         disp_b(IS) = Q(2:end,i);
                
        dis(:,i) = APDFT.Gamma*disp_b;

        vel_b(IC) = nu.*disp_b(IS);
        vel_b(IS) = -nu.*disp_b(IC);

        vel(:,i) = APDFT.Gamma*vel_b;

       
        acc_b(IC) = -nu.^2.*disp_b(IC);
        acc_b(IS) = -nu.^2.*disp_b(IS);

        acc(:,i) = APDFT.Gamma*acc_b;
        
    end
    
    fnl = Nonlinear_forces(APDFT.t, dis, vel, acc, system);                                   % produce column vectors of nonlinear forces in time-domain

    % Forward APDFT to obtain the coefficients
    Fnl = zeros(2*H+1,system.n);
    for i = 1:system.n
        Fnl(:,i) = APDFT.Gamma\fnl(:,i);
    end

    % reshape in one vector
    Fnl = reshape(transpose(Fnl),[],1);
    

end