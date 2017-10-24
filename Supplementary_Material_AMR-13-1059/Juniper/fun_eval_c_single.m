function [c1,c2] = fun_eval_c_single(param,k)

% Dispersion relation for a single shear layer between two immiscible fluids
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

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Create the coefficients of the quadratic formula
A = (1+S);
B = -2*(S*U1 + U2);
C = S*U1^2 + U2^2 - k;
B2m4AC = B.^2 - 4.*A.*C;

% Evaluate c1 and c2 with the quadratic formula
c1r = (-B + sqrt(B2m4AC))./(2*A);
c2r = (-B - sqrt(B2m4AC))./(2*A);

% Generate a stencil
stencil = logical((imag(k)>=0));

% Change the position of the branch cut
c1 = c1r.*stencil + c2r.*(1-stencil);
c2 = c2r.*stencil + c1r.*(1-stencil);

end