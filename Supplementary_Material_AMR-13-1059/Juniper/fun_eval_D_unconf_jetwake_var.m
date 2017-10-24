function [D] = fun_eval_D_unconf_jetwake_var(param,k,w)

% Dispersion relation for
% varicose perturbations of an
% unconfined jet/wake flow
% between two immiscible fluids.
% (c) Matthew Juniper 2013
%
% param is a structure containing
%   param.L   Shear ratio (U1-U2)/(U1+U2)
%   param.S   Density ratio rho1/rho2
%
% k  is a matrix contiaining wavenumbers
% c1 is a matrix containing the 1st solutions for the phase speed, omega/k
% c2 is a matrix containing the 2nd solutions for the phase speed, omega/k

% Extract the parameters
L   = param.L;
S   = param.S;
Sig = param.Sig;

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Evaluate the phase speed
c = w./k;

% Evaluate D
D = S.*(U1-c).^2.*coth(k) + (U2-c).^2 - k*Sig;

end
