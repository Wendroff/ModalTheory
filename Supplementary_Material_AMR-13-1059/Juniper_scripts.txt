function output = select_script(run, param)

% This file contains Matlab codes for local stability analysis
% using analytical dispersion relations for plug flows.
% These are parts of the tutorial:
% "Modal Stability Theory" by Matthew Juniper, Ardeshir Hanifi,
% and Vassilios Theofilis, published in Applied Mechanics Reviews,
% 66(2), 2014.

% These dispersion relations are adapted from two JFM papers:

% The Effect of Confinement on the Stability of Two-dimensional Shear flows
% M. Juniper
% Journal of Fluid Mechanics, 565, 171-195 (2006)

% The full impulse response of two-dimensional shear flows and implications for confinement
% M. Juniper
% Journal of Fluid Mechanics, 590, 163-185 (2007)
%
% (c) Matthew Juniper 2013

% The main programs are
%
% check_fun_eval_c    : Evaluate c1 and c2, given k, and check it 
%                         against the dispersion relation
% script_click.m      : Find nearest saddle point
% script_groupvel_001 : Perform a temporal stability analysis
%                         and plot the phase and group velocities
%                         of the permitted waves.
% script_spatemp_001  : Perform a spatio-temporal stability analysis
% script_spatemp_002  : Perform a spatio-temporal stability analysis
% script_spatemp_003  : Perform a spatio-temporal stability analysis
% script_spatial_001  : Perform a spatial stability anaylsis on a
%                         single shear layer
% script_spatial_002  : Perform a spatial stability anaylsis on a
%                         single shear layer both directly and via
%                         Gaster's transformation.
% script_temporal_001 : Perform a temporal stability analysis
% script_temporal_002 : Perform a temporal stability analysis for
%                         all the dispersion relations
%
% To execute, replace the argument 'run' by the string obtained
% from the corresponding function name, sans the suffix.

% fun_eval_D_* holds the dispersion relations.
% fun_eval_c_* solves the dispersion relations explicitly.
% fun_eval_D.m and fun_eval_c.m are gatekeeper functions.
% fun_eval_k.m uses an iterative method to find k(w).
% fun_eval_w.m uses an explicit method to find w(k).
% fun_eval_dwdk.m estimates the group velocity with a finite difference method.
% fun_eval_dwdk0.m uses an interative method to find k(w) at which dw/dk = 0.
% 
% Within fun_eval_c_*, the stencil changes Matlab's default ordering of the two solutions to c(k) in order to help visualization.

switch run
  case 'check_fun_eval_c'
    output = check_fun_eval_c();
  case 'script_click'
    output = script_click(param);
  case 'script_groupvel_001'
    output = script_groupvel_001();
  case 'script_spatemp_001'
    output = script_spatemp_001();
  case 'script_spatemp_002'
    output = script_spatemp_002();
  case 'script_spatemp_003'
    output = script_spatemp_003();
  case 'script_spatial_001'
    output = script_spatial_001();
  case 'script_spatial_002'
    output = script_spatial_002();
  case 'script_temporal_001'
    output = script_temporal_001();
  case 'script_temporal_002'
    output = script_temporal_002();
  otherwise
    warning('AMR:nodemo', 'No such demo');
end

end

function output = check_fun_eval_c()

% Evaluate c1 and c2, given k, and check it against the dispersion relation
% (c) Matthew Juniper 2013

clear
clc

% Input the parameters
param.L = -0.9;
param.S = 1;
param.Sig = 1.3;
param.h = 2.32;

% Input the value of k
k = complex(2*rand(1),0);

% Create a cell array containing the names of the disp. rels to check
scheme{1} = 'single';
scheme{2} = 'unconf_jetwake_sin';
scheme{3} = 'unconf_jetwake_var';
scheme{4} = 'conf_jetwake_sin';
scheme{5} = 'conf_jetwake_var';

% Cycle through the dispersion relations to check them
for nn = 1:length(scheme)
  param.scheme = scheme{nn};
  % Calculate omega1 and omega2 (w1 and w2 for short)
  [c1,c2] = fun_eval_c(param,k);
  w1 = c1*k;
  w2 = c2*k;
  % Calculate D1 and D2
  D1 = fun_eval_D(param,k,w1);
  D2 = fun_eval_D(param,k,w2);
  % Check that both are zero to machine precision
  if (abs(D1)<10*eps)&&(abs(D2)<10*eps)
    disp(['Correct dispersion relation for ',param.scheme])
  else
    disp(['Incorrect dispersion relation for ',param.scheme])
    beep
  end
end

output = param;

end

function output = script_click(param)

% With param already in memory,
% click on a point in the figure and find the nearest saddle point.
% Determine whether it is Convectively or Absolutely unstable.
% (c) Matthew Juniper 2013

% Click on a point on the figure
click = ginput(1);
k_init = complex(click(1),click(2));
tol = 10^-6;

% Converge to the nearest point, where dwdk = 0
[k,w] = fun_eval_dwdk0(param,k_init,tol);

% Plot the point accurately

% Determine whether absolutely or convectively unstable
if imag(w) > 0
  disp('Absolutely unstable')
  plot(k,'ro','MarkerFaceColor','r')
elseif imag(w) < 0
  disp('Convectively unstable')
  plot(k,'go','MarkerFaceColor','g')
else
  disp('Marginally unstable')
  plot(k,'yo','MarkerFaceColor','y')
end

output = {};

end

function output = script_groupvel_001()

% Perform a temporal stability anaylsis
% and plot the phase and group velocities of the permitted waves.
% (c) Matthew Juniper 2013

% Set the dispersion relation
param.scheme = 'single';
% Set the parameters
param.L = 1;
param.S = 1;
param.Sig = 0.001;
param.h = 1;
% Set the sampling values of the wavenumber, k
k = [0:0.01:4]';
% Evaluate the phase velocity
[c1,c2] = fun_eval_c(param,k);
% Evaluate the group velocity
[dwdk1,dwdk2] = fun_eval_dwdk(param,k);

% Plot the results on the same chart
figure(1);clf
subplot(2,2,1)
title(['Phase and group velocities '...
  ' , $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$'],'Interpreter','Latex')
%
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(k,imag(c1),'-')
plot(k,imag(c2),'-')
plot(k,imag(dwdk1),'--')
plot(k,imag(dwdk2),'--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('Imag($c_g$)','Interpreter','Latex')
set(gca,'YLim',[-4 4])
han = legend('$c_1$','$c_2$','$c_{g1}$','$c_{g2}$');
set(han,'Interpreter','Latex')
%
subplot(2,2,2)
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(k,real(c1),'-')
plot(k,real(c2),'-')
plot(k,real(dwdk1),'--')
plot(k,real(dwdk2),'--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('Real($c_g$)','Interpreter','Latex')
set(gca,'YLim',[-4 4])
han = legend('$c_1$','$c_2$','$c_{g1}$','$c_{g2}$',3);
set(han,'Interpreter','Latex')
%
% print -deps 'script_groupvel_001.eps'

output = {};

end

function output = script_spatemp_001()

% Perform a spatio-temporal stability anaylsis
% (c) Matthew Juniper 2013

param.scheme = 'single';

% Input the parameters
param.L = 0.9;
param.S = 1;
% param.Sig = 1;
% param.h = 1;

% Set the sampling values of the wavenumber, k
[kr,ki] = meshgrid(0:0.02:2,-1:0.02:1);
k = complex(kr,ki);

% Calculate omega1 and omega2 at these values of k
[w1,w2] = fun_eval_w(param,k);

% Contour plot theses side by side
figure(1);clf
% Plot the real and imaginary parts of omega1
subplot(2,2,1)
cla; hold all; set(gca,'FontSize',14,'FontName','Times')
colormap(gray)
xlabel('$k_r$','Interpreter','Latex')
ylabel('$k_i$','Interpreter','Latex')
title('$\omega_1(k)$','Interpreter','Latex')
contourf(kr,ki,imag(w1),20)
contour(kr,ki,imag(w1),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w1),20,'LineColor',[0.6 0.6 0.6])
axis image
colorbar

% Plot the real and imaginary parts of omega2
subplot(2,2,2)
cla; hold all; set(gca,'FontSize',14,'FontName','Times')
xlabel('$k_r$','Interpreter','Latex')
ylabel('$k_i$','Interpreter','Latex')
title('$\omega_2(k)$','Interpreter','Latex')
contourf(kr,ki,imag(w2),20)
contour(kr,ki,imag(w2),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w2),20,'LineColor',[0.6 0.6 0.6])
axis image
colorbar

% print -depsc2 'script_spatemp_001.eps'

output = {};

end

function output = script_spatemp_002()

% Perform a spatio-temporal stability anaylsis
% (c) Matthew Juniper 2013

clear
clc

param.scheme = 'conf_jetwake_var';

% Input the parameters
param.L = 1;
param.S = 1;
param.Sig = 1;
param.h = 1.3;

% Set the sampling values of the wavenumber, k
[kr,ki] = meshgrid(0:0.05:8,-3:0.05:1);
k = complex(kr,ki);

% Calculate omega1 and omega2 at these values of k
[c1,c2] = fun_eval_c(param,k);
w1 = c1.*k; w2 = c2.*k;

% Set the limits for w
CLim = [-1,3]; CLev = -8:0.2:8;

% Contour plot theses side by side
figure(1);clf
% Plot the real and imaginary parts of omega1
cla; hold all; set(gca,'FontSize',14,'FontName','Times')
colormap(gray)
xlabel('$k_r$','Interpreter','Latex')
ylabel('$k_i$','Interpreter','Latex')
title('$\omega_1(k)$','Interpreter','Latex')
contourf(kr,ki,imag(w1),CLev)
contour(kr,ki,imag(w1),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w1),CLev,'LineColor',[0.6 0.6 0.6])
set(gca,'CLim',CLim)

axis image

%print -dpdf 'script_spatemp_002.pdf'

output = {};

end

function output = script_spatemp_003()

% Perform a spatio-temporal stability anaylsis
% (c) Matthew Juniper 2013

clear
clc

param.scheme = 'unconf_jetwake_var';

% Input the parameters
param.L = 1;
param.S = 1;
param.Sig = 0.0;

% Set the sampling values of the wavenumber, k
[kr,ki] = meshgrid(0:0.01:2,-2:0.01:0);
k = complex(kr,ki);

% Calculate omega1 and omega2 at these values of k
[c1,c2] = fun_eval_c(param,k);
w1 = c1.*k; w2 = c2.*k;

% Set the limits for w
CLim = [-2,2]; CLev = -4:0.1:4;

% Contour plot theses side by side
figure(1);clf
% Plot the real and imaginary parts of omega1
cla; hold all; set(gca,'FontSize',14,'FontName','Times')
colormap(gray)
xlabel('$k_r$','Interpreter','Latex')
ylabel('$k_i$','Interpreter','Latex')
title('$\omega_1(k)$','Interpreter','Latex')
contourf(kr,ki,imag(w1),CLev)
contour(kr,ki,imag(w1),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w1),CLev,'LineColor',[0.6 0.6 0.6])
set(gca,'CLim',CLim)

axis image

% print -dpdf 'script_spatemp_003.pdf'

output = {};

end

function output = script_spatial_001()

% Perform a spatial stability anaylsis on a single shear layer
% (c) Matthew Juniper 2013

clear
clc

% Input the parameters
param.scheme = 'single';
param.L = 0.5;
param.S = 1;

% Set the sampling values of the angular frequency, w
w = [0.01:0.01:1]';
% Input an initial guess for k
k_init = complex(0.1,-0.02);
% Input the tolerance for convergence to the final value of k
tol = 10^-6;
% Initialize k vector
k = 0*w;

% Loop through all the values of omega, finding the relevant value of k
for nn = 1:length(w)
  k(nn) = fun_eval_k(param,k_init,w(nn),tol);
end

% Plot the results on the same chart
figure(1);clf
subplot(2,2,1)
title(['Spatial; '...
  '   $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$'],'Interpreter','Latex')
%
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),imag(k),'k')
ylabel('Imag($k$)','Interpreter','Latex')
axis([0 1 -0.15 0.05])
%
subplot(2,2,3)
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),real(k),'k')
xlabel('angular frequency, $\omega$','Interpreter','Latex')
ylabel('Real($k$)','Interpreter','Latex')
axis([0 1 0 1])

% print -deps 'script_spatial_001.eps'

output = {};

end

function output = script_spatial_002()

% Perform a spatial stability anaylsis on a single shear layer
% both directly and via Gaster's transformation.
% (c) Matthew Juniper 2013

clear
clc

% Input the parameters
param.scheme = 'single';
param.L = 0.5;
param.S = 1;

%%% TEMPORAL ANALYSIS
% Set the sampling values of the wavenumber, k
k = (0:0.01:2)';
% Evaluate the phase velocity
[w1,w2] = fun_eval_w(param,k);

%%% GASTER'S TRANSFORMATION
% Evaluate the group velocity
[dwdk1,dwdk2] = fun_eval_dwdk(param,k);
% Find kr and ki (approximately) through Gaster's transformation
krG = k;
kiG = -imag(w1)./(real(dwdk1));

%%% SPATIAL ANALYSIS (taken from script_spatial_001.m)
% Set the sampling values of the angular frequency, w
w = [0.01:0.01:1]';
% Input the tolerance for convergence to the final value of k
tol = 10^-6;
% Initialize k vector
k = 0*w;
% Loop through all the values of omega, finding the relevant value of k
for nn = 1:length(w)
  % Find the index of w1 that is closest to w
  ind = dsearchn(real(w1),w(nn));
  % Find the value of kiG at that index
  k_init = complex(krG(ind),kiG(ind));
  % Find the exact value of k from the spatial analysis
  k(nn) = fun_eval_k(param,k_init,w(nn),tol);
end

% Plot the results on the same chart
figure(1);clf
subplot(2,2,1)
title(['Spatial'...
  ' ; $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$'],'Interpreter','Latex')
%
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),imag(k),'k--')
plot(real(w1),kiG,'k')
ylabel('Imag($k$)','Interpreter','Latex')
legend('Exact','Gaster',4)
axis([0 1 -0.15 0.05])
%
subplot(2,2,3)
cla; hold on; set(gca,'FontSize',14,'FontName','Times')
plot(real(w),real(k),'k--')
plot(real(w1),krG,'k')
xlabel('angular frequency, $\omega$','Interpreter','Latex')
ylabel('Real($k$)','Interpreter','Latex')
legend('Exact','Gaster',2)
axis([0 1 0 1])

% print -deps 'script_spatial_002.eps'

output = {};

end

function output = script_temporal_001()

% Perform a temporal stability anaylsis
% (c) Matthew Juniper 2013

% Set the dispersion relation
param.scheme = 'single';
% Set the parameters
param.L = 0.9;
param.S = 1;
% Set the sampling values of the wavenumber, k
k = (0:0.01:4)';
% Evaluate the phase velocity
[w1,w2] = fun_eval_w(param,k);

% Plot the results on the same chart
figure(1);clf
subplot(2,1,1)
cla; hold on; set(gca,'FontSize',20,'FontName','Times')
plot(k,imag(w1),'k-')
plot(k,imag(w2),'k--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('growth rate, $\omega_i$','Interpreter','Latex')
title(['Single shear layer; ',...
  '$\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S),'$'],'Interpreter','Latex')
legend('First eigenvalue','Second eigenvalue')
%
subplot(2,1,2)
cla; hold on; set(gca,'FontSize',20,'FontName','Times')
plot(k,real(w1),'k-')
plot(k,real(w2),'k--')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('frequency, $\omega_r$','Interpreter','Latex')
legend('First eigenvalue','Second eigenvalue',2)
%
% print -deps 'script_temporal_001.eps'

output = {};

end

function output = script_temporal_002()

% Perform a temporal stability anaylsis for all the dispersion relations
% (c) Matthew Juniper 2013

clear
clc

% Create a cell array containing the names of the disp. rels to check
scheme{1} = 'single';
scheme{2} = 'unconf_jetwake_sin';
scheme{3} = 'unconf_jetwake_var';
scheme{4} = 'conf_jetwake_sin';
scheme{5} = 'conf_jetwake_var';

% Input the parameters
param.L = 0.9;
param.S = 1;
param.Sig = 1;
param.h = 2;

% Set the sampling values of the wavenumber, k
k = [0:0.01:2.5]';

% Set figure defaults
set(gcf,'DefaultAxesColorOrder',[0 0 0],...
  'DefaultAxesLineStyleOrder','.|o|s|+|x')

% Open the figure
figure(1);clf
subplot(1,1,1)
cla; hold all; set(gca,'FontSize',20,'FontName','Times')
xlabel('wavenumber, $k$','Interpreter','Latex')
ylabel('growth rate, $\omega_i$','Interpreter','Latex')
title([' '...
  '  $\Lambda = ',num2str(param.L), ...
  '$, $S = ',num2str(param.S), ...
  '$, $\Sigma = ',num2str(param.Sig), ...
  '$, $h = ',num2str(param.h), ...
  '$'],'Interpreter','Latex')

% Plot the symbols
for nn = 1:length(scheme)
  % Set the dispersion relation
  param.scheme = scheme{nn};
  % Evaluate the phase velocity
  [c1,c2] = fun_eval_c(param,k);
  % Evaluate the growth rates
  w1 = c1.*k; w2 = c2.*k;
  % Plot the first eigenvalue on the same chart
  plot(k(1:20:end),imag(w1(1:20:end)),'MarkerSize',6);
end

% Plot the lines
for nn = 1:length(scheme)
  % Set the dispersion relation
  param.scheme = scheme{nn};
  % Evaluate the phase velocity
  [c1,c2] = fun_eval_c(param,k);
  % Evaluate the growth rates
  w1 = c1.*k; w2 = c2.*k;
  % Plot the first eigenvalue on the same chart
  plot(k,imag(w1),'k-');
end

han = legend(scheme);
set(han,'Interpreter','none')

% print -deps2 'script_temporal_002.eps'

output = {};

end

function [c1,c2] = fun_eval_c(param,k)

% Evaluate c1 and c2, given k
% (c) Matthew Juniper 2013
%
% param.scheme contains the name of the dispersion relation

switch param.scheme
  case 'single'
    [c1,c2] = fun_eval_c_single(param,k);
  case 'unconf_jetwake_sin'
    [c1,c2] = fun_eval_c_unconf_jetwake_sin(param,k);
  case 'unconf_jetwake_var'
    [c1,c2] = fun_eval_c_unconf_jetwake_var(param,k);
  case 'conf_jetwake_sin'
    [c1,c2] = fun_eval_c_conf_jetwake_sin(param,k);
  case 'conf_jetwake_var'
    [c1,c2] = fun_eval_c_conf_jetwake_var(param,k);
  otherwise
    disp('param.scheme not recognized in fun_eval_c')
end % end switch

end % end function

function [c1,c2] = fun_eval_c_conf_jetwake_sin(param,k)

% Dispersion relation for
% sinuous perturbations of an
% confined jet/wake flow
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
h   = param.h;

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Create the coefficients of the quadratic formula
A = (coth(k*h)+S*tanh(k));
B = -2*(S*tanh(k)*U1 + coth(k*h)*U2);
C = S*U1^2*tanh(k) + U2^2*coth(k*h) - k*Sig;

B2m4AC = B.^2 - 4.*A.*C;

% Evaluate c1 and c2 with the quadratic formula
c1r = (-B + sqrt(B2m4AC))./(2*A);
c2r = (-B - sqrt(B2m4AC))./(2*A);

% Generate a stencil
stencil = logical(imag(B2m4AC)>=0);

% Change the position of the branch cut
c1 = c1r.*stencil + c2r.*(1-stencil);
c2 = c2r.*stencil + c1r.*(1-stencil);

end

function [c1,c2] = fun_eval_c_conf_jetwake_var(param,k)

% Dispersion relation for
% varicose perturbations of an
% confined jet/wake flow
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
h   = param.h;

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Create the coefficients of the quadratic formula
A = (coth(k*h)+S*coth(k));
B = -2*(S*coth(k)*U1 +coth(k*h)*U2);
C = S*U1^2*coth(k) + U2^2*coth(k*h) - k*Sig;

B2m4AC = B.^2 - 4.*A.*C;

% Evaluate c1 and c2 with the quadratic formula
c1r = (-B + sqrt(B2m4AC))./(2*A);
c2r = (-B - sqrt(B2m4AC))./(2*A);

% Generate a stencil
stencil = logical(imag(B2m4AC)>=0);

% Change the position of the branch cut
c1 = c1r.*stencil + c2r.*(1-stencil);
c2 = c2r.*stencil + c1r.*(1-stencil);

end

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

function [c1,c2] = fun_eval_c_unconf_jetwake_sin(param,k)

% Dispersion relation for
% sinuous perturbations of an
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

% Create the coefficients of the quadratic formula
A = (1+S*tanh(k));
B = -2*(S*tanh(k)*U1 + U2);
C = S*U1^2*tanh(k) + U2^2 - k*Sig;

B2m4AC = B.^2 - 4.*A.*C;

% Evaluate c1 and c2 with the quadratic formula
c1r = (-B + sqrt(B2m4AC))./(2*A);
c2r = (-B - sqrt(B2m4AC))./(2*A);

% Generate a stencil
stencil = logical(imag(B2m4AC)>=0);

% Change the position of the branch cut
c1 = c1r.*stencil + c2r.*(1-stencil);
c2 = c2r.*stencil + c1r.*(1-stencil);

end

function [c1,c2] = fun_eval_c_unconf_jetwake_var(param,k)

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

% Create the coefficients of the quadratic formula
A = (1+S*coth(k));
B = -2*(S*coth(k)*U1 + U2);
C = S*U1^2*coth(k) + U2^2 - k*Sig;

B2m4AC = B.^2 - 4.*A.*C;

% Evaluate c1 and c2 with the quadratic formula
c1r = (-B + sqrt(B2m4AC))./(2*A);
c2r = (-B - sqrt(B2m4AC))./(2*A);

% Generate a stencil
stencil = logical(imag(B2m4AC)>=0);

% Change the position of the branch cut
c1 = c1r.*stencil + c2r.*(1-stencil);
c2 = c2r.*stencil + c1r.*(1-stencil);

end

function [D] = fun_eval_D(param,k,w)

% Evaluate D given w and k
% (c) Matthew Juniper 2013
%
% param.scheme contains the name of the dispersion relation

switch param.scheme
  case 'single'
    [D] = fun_eval_D_single(param,k,w);
  case 'unconf_jetwake_sin'
    [D] = fun_eval_D_unconf_jetwake_sin(param,k,w);
  case 'unconf_jetwake_var'
    [D] = fun_eval_D_unconf_jetwake_var(param,k,w);
  case 'conf_jetwake_sin'
    [D] = fun_eval_D_conf_jetwake_sin(param,k,w);
  case 'conf_jetwake_var'
    [D] = fun_eval_D_conf_jetwake_var(param,k,w);
  otherwise
    disp('param.scheme not recognized in fun_eval_D')
end % end switch

end % end function

function [D] = fun_eval_D_conf_jetwake_sin(param,k,w)

% Dispersion relation for
% sinuous perturbations of an
% confined jet/wake flow
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
h   = param.h;

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Evaluate the phase speed
c = w./k;

% Evaluate D
D = S.*(U1-c).^2.*tanh(k) + (U2-c).^2.*coth(k*h) - k*Sig;

end

function [D] = fun_eval_D_conf_jetwake_var(param,k,w)

% Dispersion relation for
% varicose perturbations of an
% confined jet/wake flow
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
h   = param.h;

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Evaluate the phase speed
c = w./k;

% Evaluate D
D = S.*(U1-c).^2.*coth(k) + (U2-c).^2.*coth(k*h) - k*Sig;

end

function [D] = fun_eval_D_single(param,k,w)

% Dispersion relation for a single shear layer between two immiscible fluids
% (c) Matthew Juniper 2013
%
% param is a structure containing
%   param.L   Shear ratio (U1-U2)/(U1+U2)
%   param.S   Density ratio rho1/rho2
%
% k  is a matrix containing wavenumbers
% c1 is a matrix containing the 1st solutions for the phase speed, omega/k
% c2 is a matrix containing the 2nd solutions for the phase speed, omega/k

% Extract the parameters
L   = param.L;
S   = param.S;

% Create the non-dimensional velocities
U1 = 1+L;
U2 = 1-L;

% Evaluate the phase speed
c = w./k;

% Evaluate D
D = S.*(U1-c).^2 + (U2-c).^2 - k;

end

function [D] = fun_eval_D_unconf_jetwake_sin(param,k,w)

% Dispersion relation for
% sinuous perturbations of an
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
D = S.*(U1-c).^2.*tanh(k) + (U2-c).^2 - k*Sig;

end

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

function [dwdk1,dwdk2] = fun_eval_dwdk(param,k)

% Evaluate dw1/dk and dw2/dk by finite difference, given k
% (c) Matthew Juniper 2013
%
% (n.b. this version will break if k is on a branch cut)

% Input incremental k
delta = 1e-6;

% Calculate w1 and w2 at these values of k
[w11,w21] = fun_eval_w(param,k);
[w12,w22] = fun_eval_w(param,k+delta);

% Calculate the gradients
dwdk1 = (w12-w11)/delta;
dwdk2 = (w22-w21)/delta;

end % end function

function [k,w] = fun_eval_dwdk0(param,k_init,tol)

% Find the value of k that gives dw/dk = 0, given k_init
% (c) Matthew Juniper 2013

% Set the options for fsolve
options = optimset('Display','Off','TolFun',tol);

% Converge to the exact value of k using Matlab's fsolve
k = fsolve(@myfun,k_init,options);

% Evaluate omega at this point
[dwdk01,dwdk02] = fun_eval_dwdk(param,k);
[w01,w02] = fun_eval_w(param,k);
[~,i] = min(abs([dwdk01,dwdk02]));
w = (i==1)*w01 + (i==2)*w02;

% Embedded function to evaluate dwdk
  function [dwdk] = myfun(k)
    [dwdk1,dwdk2] = fun_eval_dwdk(param,k);
    dwdk = min(dwdk1,dwdk2);
  end
% End of embedded function

end % end function

function [k_out] = fun_eval_k(param,k_init,w,tol)

% Iterate towards k, given omega
% (c) Matthew Juniper 2013
%
% param.scheme is the name of the dispersion relation
% k_init       is the initial guess for k (complex)
% tol          is the tolerance required for the iteration

% Set the options for fsolve
options = optimset('Display','Off','TolFun',tol);

% Converge to the exact value of k using Matlab's fsolve
[k_out] = fsolve(@fun_eval_D_loc,k_init,options);

% Embedded function to call the correct dispersion relation
  function [D] = fun_eval_D_loc(k)
    switch param.scheme
      case 'single'
        [D] = fun_eval_D_single(param,k,w);
      case 'unconf_jetwake_sin'
        [D] = fun_eval_D_unconf_jetwake_sin(param,k,w);
      case 'unconf_jetwake_var'
        [D] = fun_eval_D_unconf_jetwake_var(param,k,w);
      case 'conf_jetwake_sin'
        [D] = fun_eval_D_conf_jetwake_sin(param,k,w);
      case 'conf_jetwake_var'
        [D] = fun_eval_D_conf_jetwake_var(param,k,w);
      otherwise
        disp('param.scheme not recognized in fun_eval_D')
    end % end switch
  end
% End of embedded function

end % end function

function [w1,w2] = fun_eval_w(param,k)

% Evaluate w1 and w2, given k
% (c) Matthew Juniper 2013
%
% param.scheme contains the name of the dispersion relation

switch param.scheme
  case 'single'
    [c1,c2] = fun_eval_c_single(param,k);
  case 'unconf_jetwake_sin'
    [c1,c2] = fun_eval_c_unconf_jetwake_sin(param,k);
  case 'unconf_jetwake_var'
    [c1,c2] = fun_eval_c_unconf_jetwake_var(param,k);
  case 'conf_jetwake_sin'
    [c1,c2] = fun_eval_c_conf_jetwake_sin(param,k);
  case 'conf_jetwake_var'
    [c1,c2] = fun_eval_c_conf_jetwake_var(param,k);
  otherwise
    disp('param.scheme not recognized in fun_eval_c')
end % end switch

% Evaluate w1 and w2 from c1 and c2
w1 = k.*c1; w2 = k.*c2;

end % end function