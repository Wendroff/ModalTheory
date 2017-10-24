function output = script_spatemp_003()

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
[kr,ki] = meshgrid(-1:0.05:3,-7:0.05:1);
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
% contourf(kr,ki,imag(w1),CLev)
contour(kr,ki,imag(w1),[0 0],'LineColor','k','LineWidth',2)
contour(kr,ki,real(w1),CLev,'LineColor',[0.6 0.6 0.6])
set(gca,'CLim',CLim)

axis image

% print -dpdf 'script_spatemp_003.pdf'
figure;
surf(kr,ki,imag(w1))
axis([-1 3 -7 1 -100 100])
output = {};

end
